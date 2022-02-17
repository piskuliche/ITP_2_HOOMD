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

def set_snap(sselem,ff,ffelem,nmols):
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
    sselem.N = len(ffelem)*nmols
    num_atoms = len(ff.amap)
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

    # make copies for each molecule
    natms = len(ff.amap)
    el = np.array(list(map(type_idlist,elemtypes)))
    og = np.array(numgroups)
    outelem, outgroup = el,og
    for i in range(nmols-1):
        outelem = np.concatenate((outelem,el),axis=0)
        outgroup = np.concatenate((outgroup,np.add(og,natms*(i+1))),axis=0)
    #exit()
    sselem.typeid = outelem
    sselem.group = outgroup
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

def Generate_Snapshot(mol_ff,forcefield,sysinfo):
    """
    This function creates an initial configuration and builds a gsd snapshot object. It writes it to the file 'molecular.gsd', which can then be saved.
    Inputs:
        -mol_ff: molecule_info class object
        -forcefield: forcefield_info class object
        -sysinfo: dictionary with molecule names as keys, and numbers of molecules.
    """
    # Converts charges to proper internal units
    charge_conv = 0.0848385920

    # Initialize a new snapshot
    snapshot = gsd.hoomd.Snapshot()

    # Set number of atoms based off of gro file
    snapshot.particles.N = N

    # Initialize particle positions
    snapshot.particles.position = r
    setcharge = lambda i: mol_ff.atoms[i].charge*charge_conv
    snapshot.particles.charge = np.array(list(map(setcharge,names)))

    # Initialize particle types, find unique types, and assign type ids
    # These use the map function to get arrays of all the types
    snapshot.particles = set_part(snapshot.particles,ff["DOPC"])

    # Set Box Size
    snapshot.configuration.box = [200, 200, 200, 0, 0, 0]

    # Set bonds,angles,dihedrals between atoms
    snapshot.pairs = set_snap(snapshot.pairs,mol_ff, mol_ff.bonds,sysinfo["DOPC"])
    snapshot.bonds = set_snap(snapshot.bonds, mol_ff, mol_ff.bonds,sysinfo["DOPC"])
    snapshot.angles  = set_snap(snapshot.angles, mol_ff, mol_ff.angles,sysinfo["DOPC"])
    snapshot.dihedrals = set_snap(snapshot.dihedrals, mol_ff, mol_ff.dihedrals,sysinfo["DOPC"])

    # Save snapshot
    with gsd.hoomd.open(name='molecular.gsd', mode='wb') as f:
        f.append(snapshot)

    return snapshot

def Generate_Forces(snapshot,forcefieldi,r_cut=1.2,buffer=0.1,alpha=1,coul_order=6,res=64):
    """
    This function generates the forces object by pulling together all of the forces objects and calling them from the forcefield.
    Inputs:
        -snapshot object: GSD snapshot object that has been previously generated, for example by Generate_Snapshot
        -forcefield: forcefield_info class object.
    Outputs:
        -all_forces: a list of hoomd.md.XXX objects where XXX stands for the different types of forces in the system.
    """
    nl = hoomd.md.nlist.Cell(buffer=0.1)
    lj = hoomd.md.pair.LJ(nl, default_r_cut=1.2)
    lj.params = gen_pair(forcefield)
    ew,coul = hoomd.md.long_range.pppm.make_pppm_coulomb_forces(nl, resolution=(res,res,res),order=coul_order,r_cut=r_cut)

    # Setup HOOMD Simulation
    # Special Pair Potential 
    spec_lj_forces = hoomd.md.special_pair.LJ()
    spec_coul_forces = hoomd.md.special_pair.Coulomb()
    for pair in snapshot.pairs.types:
        ptuple = tuple(pair.split("-"))
        spec_lj_forces.params[pair] = lj.params[ptuple]
        spec_lj_forces.r_cut[pair] = r_cut
        spec_coul_forces.params[pair] = dict(alpha=alpha)
        spec_coul_forces.r_cut[pair] = r_cut
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

    all_forces=[ew,coul,spec_lj_forces,spec_coul_forces,bondforces,angleforces,dihedralforces]
    return all_forces


if __name__ == "__main__":
    # Converts charges to proper internal units
    charge_conv = 0.0848385920
    # Read the # of atoms, and the positions
    N,r,names = read_gro("test.gro")
    # Read in molecular force field
    ff = {}
    ff["DOPC"] = read_mol_ff("DOPC-OPLS.itp")
    sysinfo={"DOPC":72}

    # Read in total forcefield
    forcefield = read_ff_itp("OPLSAA.itp")

    # Generate Snapshot
    snapshot=Generate_Snapshot(ff["DOPC"],forcefield,sysinfo)
   
    # Generate Forces
    all_forces = Generate_Forces(snapshot,forcefield)

    # Perform the MD simulation.
    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
    sim.create_state_from_gsd(filename='molecular.gsd')

    # Set initial velocities
    sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=2.5)
    fire = hoomd.md.minimize.FIRE(dt=0.05,
                            force_tol=1e-2,
                            angmom_tol=1e-2,
                            energy_tol=1e-7)
    fire.methods.append(hoomd.md.methods.NVE(hoomd.filter.All()))
    sim.operations.integrator = fire
    fire.converged()
    while not(fire.converged):
        sim.run(100)

    #ensemble = hoomd.md.methods.NPT(filter=hoomd.filter.All(),kT=2.5,tau=1, tauS=10, S=[1/16,1/16,1/16,0,0,0],couple='xy')
    ensemble = hoomd.md.methods.NVE(filter=hoomd.filter.All())
    integrator = hoomd.md.Integrator(dt=0.001,
                                     methods=[ensemble],
                                     forces=all_forces)
    gsd_writer = hoomd.write.GSD(filename='molecular_trajectory.gsd',
                                 trigger=hoomd.trigger.Periodic(1),
                                 mode='wb')

    # Calculate Thermodynamics
    thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
    sim.operations.computes.append(thermodynamic_properties)
    
    # Add A file to Log Thermodynamics
    logger = hoomd.logging.Logger()
    logger.add(thermodynamic_properties)

    gsd_writer2 = hoomd.write.GSD(filename='log.gsd',
                             trigger=hoomd.trigger.Periodic(1),
                             mode='wb',
                             filter=hoomd.filter.Null())

    # Final Run
    sim.operations.integrator = integrator
    sim.operations.writers.append(gsd_writer)
    sim.operations.writers.append(gsd_writer2)
    gsd_writer2.log = logger
    sim.run(10)
    print(thermodynamic_properties.kinetic_temperature)


