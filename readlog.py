import hoomd
import gsd.hoomd
import numpy as np

traj = gsd.hoomd.open('log.gsd', 'rb')
for frame in traj:
    print(frame.log['md/compute/ThermodynamicQuantities/potential_energy'])
