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

def set_snap(sselem,ffatoms,ffelem):
    sselem.N = len(ffelem)
    nset = lambda atm: list(names[np.subtract(atm,1)])
    nlist = list(map(tuple,list(map(nset,ffelem))))
    elemtypes = []
    numgroups = []
    for elem in nlist:
        tmp = []
        tmpg = []
        for atm in elem:
            tmp.append(ffatoms[atm].type)
            tmpg.append(np.where(names==atm)[0][0])
        elemtypes.append(tuple(tmp))
        numgroups.append(tmpg)
    unique, u_ids = np.unique(elemtypes,axis=0),np.arange(len(np.unique(elemtypes,axis=0)))
    utup = list(map(tuple,unique))
    type_id_map = dict(zip(utup, u_ids))
    type_idlist = lambda etype: type_id_map[etype]
    sselem.types=list(map(tuple,unique))
    sselem.typeid = np.array(list(map(type_idlist,elemtypes)))
    sselem.groups = numgroups
    return sselem


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
    typelist = lambda atm: ff["DOPC"].atoms[atm].type
    snapshot.particles.types = np.array(list(map(typelist,names)))
    unique,u_ids = np.unique(snapshot.particles.types),np.arange(len(np.unique(snapshot.particles.types)))
    type_id_map = dict(zip(unique,u_ids))
    type_idlist = lambda atype: type_id_map[atype]
    snapshot.particles.typeid = np.array(list(map(type_idlist,snapshot.particles.types)))

    # Set bonds between atoms
    """
    snapshot.bonds.N = len(ff["DOPC"].bonds)
    nset = lambda atm: list(names[np.subtract(atm,1)])
    nlist = list(map(tuple,list(map(nset,ff["DOPC"].bonds))))
    bondtypes = []
    numgroups = []
    for bond in nlist:
        tmp = []
        tmpg = []
        for atm in bond:
            tmp.append(ff["DOPC"].atoms[atm].type)
            tmpg.append(np.where(names==atm)[0][0])
        bondtypes.append(tuple(tmp))
        numgroups.append(tmpg)
    snapshot.bonds.types=bondtypes
    unique, u_ids = np.unique(snapshot.bonds.types,axis=0),np.arange(len(np.unique(snapshot.bonds.types,axis=0)))
    utup = list(map(tuple,unique))
    type_id_map = dict(zip(utup, u_ids))
    type_idlist = lambda btype: type_id_map[btype]
    snapshot.bonds.typeid = np.array(list(map(type_idlist,snapshot.bonds.types)))
    snapshot.bonds.groups = numgroups
    """
    snapshot.bonds = set_snap(snapshot.bonds, ff["DOPC"].atoms, ff["DOPC"].bonds)
    snapshot.angles  = set_snap(snapshot.angles, ff["DOPC"].atoms, ff["DOPC"].angles)
    snapshot.dihedrals = set_snap(snapshot.dihedrals, ff["DOPC"].atoms, ff["DOPC"].dihedrals)

