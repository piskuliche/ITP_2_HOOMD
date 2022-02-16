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
        -r: numpy array in order of atoms, shape=(N,3)
        -names: numpy array of atoms names, shape=(N)
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
    """
    This takes the particles gsd snapshot object, and the molecular forcefield class object,
    and sets the objects

    Inputs:
        -sspart: GSD snapshot.particles object
        -ff: molecule_info class object
    Outputs:
        -sspart: GSD snapshot.particles object
    Sets:
        -sspart.types: unique atom types present, shape=(# unique types)
        -sspart.typeid: type ids, shape=(N)
    """
    # Make a list of atom types
    typelist = lambda atm: ff.atoms[atm].type
    atomtypes = np.array(list(map(typelist,names)))
    # Fin unique atom types
    unique,u_ids = np.unique(atomtypes),np.arange(len(np.unique(atomtypes)))
    sspart.types=unique
    # Map these to integer type ids
    type_id_map = dict(zip(unique,u_ids))
    type_idlist = lambda atype: type_id_map[atype]
    sspart.typeid = np.array(list(map(type_idlist,atomtypes)))
    return sspart

def set_snap(sselem,ff,ffelem):
    """
    This code takes the molecular info and the forcefield info and sets the bonds/angles/dihedrals
    Input:
        -sselem: snapshot.[bonds,angles, or dihedrals] object
        -ff: molecule_info class object (e.g. ff["DOPC"])
        -ffelem: corresponding piece (e.g. ff["DOPC"].bonds)
    Output:
        -sselem: Returns set snapshot.[bonds,angles,dihedrals] object
    Sets:
        -sselem.N: Number of these objects in ffelem. 
        -sselem.types: Numpy array of unique [bond,angle,dihedral] types found, shape=(# unique)
        -sselem.typeid: Numpy array of corresponding integer ids to unique types, shape=(# unique)
        -sselem.group: Numpy array of integer atom ids for each elem, shape=(sselem.N,#atomsperelem)
            Note: #atomsperelem is 2 for bonds, 3 for angles, 4 for dihedrals
    ToDo:
        -Need to modify sselem.N to account for multiple copies of each molecule
        -Need to do the same for sselem.group
        -tmpg also probably breaks for more than 1 molecule, due to names being the same between molecules
    """
    #Set # of bonds,angles, or dihedrals
    sselem.N = len(ffelem)
    elemtypes,numgroups = [],[]
    # loop over bonds
    for elem in ffelem:
        tmp,tmpg = "",[]
        # loop over atoms in each bond
        for atom in elem:
            # Get atom name
            atm = ff.amap[atom]
            # Append to string and include "-"
            tmp += ff.atoms[atm].type+"-"
            # Append location where names == atm to tmpg
            tmpg.append(np.where(names==atm)[0][0])
            print(np.where(names==atm))
        # Remove last "-"
        tmp = tmp[:-1]
        # Append string name to elemtypes
        elemtypes.append(tmp)
        # Append location to numgroups
        numgroups.append(tmpg)
    # Find Unique Types
    unique, u_ids = np.unique(elemtypes),np.arange(len(np.unique(elemtypes)))
    type_id_map = dict(zip(unique, u_ids))
    type_idlist = lambda etype: type_id_map[etype]
    sselem.types=np.array(unique)
    sselem.typeid = np.array(list(map(type_idlist,elemtypes)))
    # Set group
    sselem.group = np.array(numgroups)
    return sselem

def RB_to_OPLS(df):
    """
    This function takes a dihedraltype class object (from the forcefield_info class object)
    and converts it from RB type dihedrals to a more traditional fourier dihedral.

    This should be moved into the forcefield_info dihedraltype class as a class function. 
    """
    c0,c1,c2,c3,c4,c6=df.C
    F1 = 2*c2 + 2*c0 -2*c4-c3/2
    F2 = c4 - c2 + c3/2
    F3 = -c3/2
    F4 = -c4/4
    return dict(k1=F1,k2=F2,k3=F3,k4=F4)

def gen_pair(ff):
    """
    This function takes a forcefield_info class object and pulls the lj parameters from it. It also applies the geometric mixing rules to all parameters.

    Inputs: 
        -ff: forcefield_info class object
    Outputs:
        -params: Dictionary (keyed by tuple of unique atom combos) of dictionaries (two keys: sigma, epsilon)

    Future:
        Could implement options for mixing rules (for instance aritmetic)
        Should probably make it clearer between functions what ff is
    """
    def mix_geometric(xi,xj):
        # Applies geometric mixing rules.
        return np.sqrt(xi*xj)

    # Generate list of unique atom types
    atoms = []
    for key in ff.sigma:
        atoms.append(key)
    atoms = np.unique(atoms)
    params={}
    # Loop over atoms
    for ai in atoms:
        # Loop over atoms
        for aj in atoms:
            # e.g. entry=(N,O)
            #params[(N,O)] = {"sigma":sigma_{NO}, "epsilon":epsilon_{NO}}
            entry = tuple([ai,aj])
            params[entry]={"sigma":mix_geometric(ff.sigma[ai],ff.sigma[aj]),"epsilon":mix_geometric(ff.epsilon[ai],ff.sigma[aj])}
    return params
                

if __name__ == "__main__":
    # Converts charges to proper internal units
    charge_conv = 0.0848385920

    ff ={}

    # Initialize a new snapshot
    snapshot = gsd.hoomd.Snapshot()
    
    # Read the # of atoms, and the positions
    N,r,names = read_gro("test.gro")

    # Read in molecular force field
    ff["DOPC"] = read_mol_ff("DOPC-OPLS.itp")

    # Read in total forcefield
    forcefield = read_ff_itp("OPLSAA.itp")

    # Set number of atoms based off of gro file
    snapshot.particles.N = N

    # Initialize particle positions
    snapshot.particles.position = r
    setcharge = lambda i: ff["DOPC"].atoms[i].charge*charge_conv
    snapshot.particles.charge = np.array(list(map(setcharge,names)))
    
    # Initialize particle types, find unique types, and assign type ids
    # These use the map function to get arrays of all the types
    snapshot.particles = set_part(snapshot.particles,ff["DOPC"]) 

    # Set Box Size
    snapshot.configuration.box = [100, 100, 100, 0, 0, 0]

    # Set bonds,angles,dihedrals between atoms
    snapshot.bonds = set_snap(snapshot.bonds, ff["DOPC"], ff["DOPC"].bonds)
    snapshot.angles  = set_snap(snapshot.angles, ff["DOPC"], ff["DOPC"].angles)
    snapshot.dihedrals = set_snap(snapshot.dihedrals, ff["DOPC"], ff["DOPC"].dihedrals)
    # Save snapshot
    with gsd.hoomd.open(name='molecular.gsd', mode='wb') as f:
        f.append(snapshot)

    # Below this point, only FF information should be generated. 
    nl = hoomd.md.nlist.Cell(buffer=0.4)
    lj = hoomd.md.pair.LJ(nl, default_r_cut=1.2)
    lj.params = gen_pair(forcefield)
    ew,coul = hoomd.md.long_range.pppm.make_pppm_coulomb_forces(nl, resolution=(64,64,64),order=6,r_cut=1.2)

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
    ensemble = hoomd.md.methods.NVE(filter=hoomd.filter.All())
    integrator = hoomd.md.Integrator(dt=0.001,
                                     methods=[ensemble],
                                     forces=[lj,ew,coul,bondforces,angleforces,dihedralforces])
    gsd_writer = hoomd.write.GSD(filename='molecular_trajectory.gsd',
                                 trigger=hoomd.trigger.Periodic(1),
                                 mode='xb')
    sim.operations.integrator = integrator
    sim.operations.writers.append(gsd_writer)
    sim.run(10)


