import hoomd
import gsd.hoomd
import numpy as np
from itp_to_hoomd import read_ff_itp,read_mol_ff,molecule_info,forcefield_info

def read_gro(fname):
    """
    This function takes a file name for a gro file, and reads it.

    Inputs: 
        -filename: filename to be read
    Outputs: 
        -N: integer number of atoms
        -r: numpy array in order of atoms
        -names: numpy array of atoms names
    """
    f=open(fname,'r')
    f.readline()
    r = []

    N=int(f.readline().strip())
    names=[]
    for i in range(N):
        _, name, _, x, y, z = f.readline().strip().split()
        r.append([float(x),float(y),float(z)])
        names.append(name)
    return N,np.array(r),np.array(names)

def set_part(sspart,ff):
    typelist = lambda atm: ff.atoms[atm].type
    atomtypes = np.array(list(map(typelist,names)))
    unique,u_ids = np.unique(atomtypes),np.arange(len(np.unique(atomtypes)))
    sspart.types=unique
    type_id_map = dict(zip(unique,u_ids))
    type_idlist = lambda atype: type_id_map[atype]
    sspart.typeid = np.array(list(map(type_idlist,atomtypes)))
    return sspart

def set_snap(sselem,ff,ffelem):
    sselem.N = len(ffelem)
    elemtypes = []
    numgroups = []
    for elem in ffelem:
        tmp = ""
        tmpg = []
        for atom in elem:
            atm = ff.amap[atom]
            tmp += ff.atoms[atm].type+"-"
            tmpg.append(np.where(names==atm)[0][0])
        tmp = tmp[:-1]
        elemtypes.append(tmp)
        numgroups.append(tmpg)
    unique, u_ids = np.unique(elemtypes),np.arange(len(np.unique(elemtypes)))
    type_id_map = dict(zip(unique, u_ids))
    type_idlist = lambda etype: type_id_map[etype]
    sselem.types=np.array(unique)
    sselem.typeid = np.array(list(map(type_idlist,elemtypes)))
    sselem.group = np.array(numgroups)
    return sselem

def RB_to_OPLS(df):
    c0,c1,c2,c3,c4,c6=df.C
    F1 = 2*c2 + 2*c0 -2*c4-c3/2
    F2 = c4 - c2 + c3/2
    F3 = -c3/2
    F4 = -c4/4
    return dict(k1=F1,k2=F2,k3=F3,k4=F4)

def gen_pair(ff):
    def mix_geometric(xi,xj):
        return np.sqrt(xi*xj)
    atoms = []
    for key in ff.sigma:
        atoms.append(key)
    atoms = np.unique(atoms)
    print(atoms)
    params={}
    for ai in atoms:
        for aj in atoms:
            entry = tuple([ai,aj])
            params[entry]={"sigma":mix_geometric(ff.sigma[ai],ff.sigma[aj]),"epsilon":mix_geometric(ff.epsilon[ai],ff.sigma[aj])}
    return params
                


if __name__ == "__main__":
    ff ={}

    # Initialize a new snapshot
    snapshot = gsd.hoomd.Snapshot()
    
    # Read the # of atoms, and the positions
    N,r,names = read_gro("DOPC-OPLS.gro")

    # Read in molecular force field
    ff["DOPC"] = read_mol_ff("DOPC-OPLS.itp")

    # Read in total forcefield
    forcefield = read_ff_itp("OPLSAA.itp")

    # Set number of atoms based off of gro file
    snapshot.particles.N = N

    # Initialize particle positions
    snapshot.particles.position = r
    
    # Initialize particle types, find unique types, and assign type ids
    # These use the map function to get arrays of all the types
    snapshot.particles = set_part(snapshot.particles,ff["DOPC"]) 
    print(snapshot.particles.types)

    # Set Box Size
    snapshot.configuration.box = [100, 100, 100, 0, 0, 0]

    # Set bonds,angles,dihedrals between atoms
    snapshot.bonds = set_snap(snapshot.bonds, ff["DOPC"], ff["DOPC"].bonds)
    snapshot.angles  = set_snap(snapshot.angles, ff["DOPC"], ff["DOPC"].angles)
    snapshot.dihedrals = set_snap(snapshot.dihedrals, ff["DOPC"], ff["DOPC"].dihedrals)
    # Save snapshot
    with gsd.hoomd.open(name='molecular.gsd', mode='wb') as f:
        f.append(snapshot)
    nl = hoomd.md.nlist.Cell(buffer=0.4)
    lj = hoomd.md.pair.LJ(nl, default_r_cut=1.2)
    lj.params = gen_pair(forcefield)

    # Setup HOOMD Simulation
    # Bond Potential
    bondforces = hoomd.md.bond.Harmonic()
    for bond in snapshot.bonds.types:
        btuple = tuple(bond.split("-"))
        bvals = forcefield.bondtypes[btuple]
        bondforces.params[bond] = dict(k=bvals.K, r0=bvals.r0)
    angleforces = hoomd.md.angle.Harmonic()
    # Angle Potential
    for angle in snapshot.angles.types:
        atuple = tuple(angle.split("-"))
        avals = forcefield.angletypes[atuple]
        angleforces.params[angle] = dict(k=avals.Kth, t0=np.radians(avals.th0))
    dihedralforces = hoomd.md.dihedral.OPLS()
    # Dihedral Potential
    for dih in snapshot.dihedrals.types:
        dtuple = tuple(dih.split("-"))
        dvals = forcefield.dihedraltypes[dtuple]
        dihedralforces.params[dih] = RB_to_OPLS(dvals)

    # Perform the MD simulation.
    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
    sim.create_state_from_gsd(filename='molecular.gsd')
    langevin = hoomd.md.methods.NVE(filter=hoomd.filter.All())
    integrator = hoomd.md.Integrator(dt=0.001,
                                     methods=[langevin],
                                     forces=[lj,bondforces,angleforces,dihedralforces])
    gsd_writer = hoomd.write.GSD(filename='molecular_trajectory.gsd',
                                 trigger=hoomd.trigger.Periodic(10),
                                 mode='xb')
    sim.operations.integrator = integrator
    sim.operations.writers.append(gsd_writer)
    sim.run(500)


