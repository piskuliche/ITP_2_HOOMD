import numpy as np

h=open("test.gro",'w')
x,y,z=[],[],[]
with open("DOPC-OPLS.gro",'r') as f:
    with open("step7_1.gro",'r') as g:
        f.readline()
        N=int(f.readline().strip())
        g.readline()
        g.readline()
        h.write("GROMACS\n")
        h.write("%d\n" % N)
        for i in range(N):
            line = f.readline().strip().split()
            gline = g.readline().strip().split()
            x.append(float(gline[3]))
            y.append(float(gline[4]))
            z.append(float(gline[5]))
            gline[1]=line[1]
            h.write("%s %s %d %10.5f %10.5f %10.5f\n" % (gline[0],gline[1],int(gline[2]),float(gline[3]),float(gline[4]),float(gline[5])))
        h.write("0.0 0.0 0.0\n")
    print(np.min(x),np.min(y),np.min(z))
    print(np.max(x),np.max(y),np.max(z))
