import numpy as np

h=open("test.gro",'w')
N=0
n=72
names=[]
with open("DOPC-OPLS.gro",'r') as f:
    f.readline()
    N=int(f.readline().strip())
    for i in range(N):
        _,name,_,_,_,_ = f.readline().strip().split()
        name = str(name)
        names.append(name)
xs,ys,zs=[],[],[]
with open("step5_input.gro",'r') as g:
    g.readline()
    g.readline()
    h.write("GROMAC\n")
    h.write("%d\n" %(N*n))
    for i in range(n):
        for j in range(N):
            m, n, aid, x, y, z = g.readline().strip().split()
            #m, n, aid, x, y, z,vx,vy,vz = g.readline().strip().split()
            xs.append(float(x))
            ys.append(float(y))
            zs.append(float(z))
            h.write("%s %s %d %10.5f %10.5f %10.5f\n" % (m,names[j],int(aid),float(x),float(y),float(z)))
print(np.max(xs)-np.min(xs))
print(np.max(ys)-np.min(ys))
print(np.max(zs)-np.min(zs))
