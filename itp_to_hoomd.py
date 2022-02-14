import numpy as np
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('-ff', default="forcefield.itp",type=str,help='Forcefield file')
parser.add_argument('-mf', default="DOPC.itp",type=str,help='Molecule file')
args = parser.parse_args()

mfile = args.mf
ffile = args.ff

class molecule_info:
    """
    This class stores the info from a molecular itp file.
    """
    class atom_info:
        def __init__(self,Q,M,TY,aid):
            self.charge = float(Q)
            self.mass = float(M)
            self.type = TY
            self.id = int(aid)
    def __init__(self):
        #atoms
        self.atoms = {}
        self.amap = {}
        #bonds
        self.bonds = []
        self.bondfunc = []
        #pairs
        self.pairs = []
        self.pairfunc = []
        # angles
        self.angles = []
        self.anglefunc = []
        # dihedrals
        self.dihedrals = []
        self.dihfunc = []
        return
    def read_moleculetype(self,line):
        self.name,self.nrexcl = line
        return
    def read_atoms(self,line):
        aid,types,_,_,names,_,charges,masses = line
        self.atoms[names]=self.atom_info(charges,masses,types,aid)
        self.amap[int(aid)] = names
        return
    def read_bonds(self,line):
        ai, aj, bondfunc = line
        ai,aj = int(ai),int(aj)
        self.bonds.append([ai,aj])
        self.bondfunc.append(bondfunc)
        return
    def read_pairs(self,line):
        ai, aj, pairfunc = line
        ai,aj = int(ai),int(aj)
        self.pairs.append([ai,aj])
        self.bondfunc.append(pairfunc)
        return
    def read_angles(self,line):
        ai, aj, ak, anglefunc = line
        ai,aj,ak = int(ai),int(aj),int(ak)
        self.angles.append([ai,aj,ak])
        self.anglefunc.append(anglefunc)
        return
    def read_dihedrals(self,line):
        ai,aj,ak,al,dihfunc = line
        ai,aj,ak,al = int(ai),int(aj),int(ak),int(al)
        self.dihedrals.append([ai,aj,ak,al])
        self.dihfunc.append(dihfunc)
        return

class forcefield_info:
    class bondtype:
        def __init__(self, r0,K):
            self.K = float(K)
            self.r0 = float(r0)
            return
    class pairtype:
        def __init__(self,sigma,epsilon):
            self.sigma = float(sigma)
            self.epsilon = float(epsilon)
            return
    class angletype:
        def __init__(self,th0,Kth):
            self.th0 = float(th0)
            self.Kth = float(Kth)
            return
    class dihedraltype:
        def __init__(self,C):
            self.C = list(map(float,C))
            return
    def __init__(self):
        self.sigma = {}
        self.epsilon = {}
        self.nbpair, self.nbsigma, self.nbepsilon = [],[],[]
        self.bondtypes = {}
        self.pairtypes = {}
        self.angletypes = {}
        self.dihedraltypes = {}
        return
    def add_lj(self,line):
        name, _,_,_,sigma,epsilon = line
        self.sigma[name] = float(sigma)
        self.epsilon[name] = float(epsilon)
        return
    def add_nonbond(self,line):
        ai, aj, _, sigma, epsilon = line
        self.nbpair.append([ai,aj])
        self.nbsigma.append(sigma)
        self.nbepsilon.append(epsilon)
        return
    def add_bondtypes(self,line):
        ai,aj,_,r0,K = line
        self.bondtypes[(ai,aj)]=self.bondtype(r0,K)
        self.bondtypes[(aj,ai)]=self.bondtype(r0,K)
        return
    def add_pairtypes(self,line):
        ai,aj,_,sigma,epsilon=line
        self.pairtypes[(ai,aj)]=self.pairtype(sigma,epsilon)
        self.pairtypes[(aj,ai)]=self.pairtype(sigma,epsilon)
        return
    def add_angletypes(self,line):
        ai,aj,ak,_,th0,Kth = line
        self.angletypes[(ai,aj,ak)]=self.angletype(th0,Kth)
        self.angletypes[(ak,aj,ai)]=self.angletype(th0,Kth)
        return
    def add_dihedraltypes(self,line):
        ai,aj,ak,al,_,C0,C1,C2,C3,C4,C5 = line
        self.dihedraltypes[(ai,aj,ak,al)]=self.dihedraltype([C0,C1,C2,C3,C4,C5])
        self.dihedraltypes[(al,ak,aj,ai)]=self.dihedraltype([C0,C1,C2,C3,C4,C5])
        return




def read_mol_ff(fname):
    """
    Reads the molecule itp file
    """
    f=open(fname,'r')
    molecule = molecule_info()
    lines = f.readlines()
    flag=""
    for line in lines:
        if ";" in line: line = line.split(";")[0]
        if "[" in line:
            flag = line.strip().split()[1]
            continue
        line = line.strip().split()
        if flag == "moleculetype" and len(line) == 2:
            molecule.read_moleculetype(line)
        if flag == "atoms" and len(line) == 8:
            molecule.read_atoms(line)
        if flag == "bonds" and len(line) == 3:
            molecule.read_bonds(line)
        if flag == "pairs" and len(line) == 3:
            molecule.read_pairs(line)
        if flag == "angles" and len(line) == 4:
            molecule.read_angles(line)
        if flag == "dihedrals" and len(line) == 5:
            molecule.read_dihedrals(line)

    print("There are %d atoms" % len(molecule.atoms))
    print("There are %d bonds" % len(molecule.bonds))
    print("There are %d pairs" % len(molecule.pairs))
    print("There are %d angles" % len(molecule.angles))
    print("There are %d dihedrals" % len(molecule.dihedrals))
    f.close()
    return molecule

def read_ff_itp(fname):
    """
    Reads the force field file
    """
    f=open(fname,'r') 
    ff = forcefield_info()
    lines = f.readlines()
    flag=""
    for line in lines:
        if ";" in line: line = line.split(";")[0]
        if "[" in line:
            flag = line.strip().split()[1]
        line = line.strip().split()
        if flag == "atomtypes" and len(line) == 6:
            ff.add_lj(line)
        if flag == "nonbond_params" and len(line) == 5:
            ff.add_nonbond(line)
        if flag == "bondtypes" and len(line) == 5:
            ff.add_bondtypes(line)
        if flag == "pairtypes" and len(line) == 5:
            ff.add_pairtypes(line)
        if flag == "angletypes" and len(line) == 6:
            ff.add_angletypes(line)
        if flag == "dihedraltypes" and len(line) == 11:
            ff.add_dihedraltypes(line)
    print("There are %d defined atoms" % len(ff.sigma))
    print("There are %d bond types" % len(ff.bondtypes))
    print("There are %d angle types" % len(ff.angletypes))
    print("There are %d dihedral types" % len(ff.dihedraltypes))
    f.close()
    return ff



if __name__ == "__main__":
    DOPC = read_mol_ff(mfile)
    ff       = read_ff_itp(ffile)
