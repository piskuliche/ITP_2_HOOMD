import gsd.hoomd,hoomd
import numpy as np

frames = gsd.hoomd.open("molecular_trajectory.gsd",'rb')

with gsd.hoomd.open(name="test.gsd",mode='wb') as g:
    for frame in frames:
        tmpframe = frame
        tmpframe.particles.position*=10
        g.append(tmpframe)

