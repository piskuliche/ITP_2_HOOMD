import numpy as np

h=open("test.gro",'w')
N=0
n=2
names=[]
with open("DOPC-OPLS.gro",'r') as f:
    f.readline()
    N=int(f.readline().strip())
    for i in range(N):
        _,name,_,_,_,_ = f.readline().strip().split()
        name = str(name)
        names.append(name)

with open("step7_1.gro",'r') as g:
    g.readline()
    g.readline()
    h.write("GROMAC\n")
    h.write("%d\n" %(N*n))
    for i in range(n):
        for j in range(N):
            m, n, aid, x, y, z,vx,vy,vz = g.readline().strip().split()
            h.write("%s %s %d %10.5f %10.5f %10.5f\n" % (m,names[j],int(aid),float(x),float(y),float(z)))
        
