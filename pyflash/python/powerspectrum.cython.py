#!./py.sh

import numpy as np
import sys

from pyflash import FlashFile
import vectoranalysis3D as va3D
from shell_avg import shell_avg_3d

def product(it):
    p = 1
    for i in it:
        p *= i
    return p

try:
    fp = sys.argv[1]
except IndexError:
    print("usage: %s <checkpoint file>" % sys.argv[0])
    sys.exit(1)

flash = FlashFile(fp)

time = flash.meta['real scalars']['time']
step = flash.meta['integer scalars']['nstep']

c_s  = flash.meta['real runtime']['c_ambient']
rho0 = flash.meta['real runtime']['rho_ambient']

GS = flash.meta['grid size']
CS = flash.meta['cell size']
DS = flash.meta['domain size']

Vgrid   = product(GS)
Vcell   = product(CS) 
Vdomain = product(DS) 

density   = flash.get_box('dens')
velocity  = [flash.get_box('vel'+dim) for dim in 'x y z'.split()]
#vorticity = va3D.curl(velocity[0],velocity[1],velocity[2],CS[0],CS[1],CS[2])

ekin = density/2.0 * va3D.norm(velocity[0],velocity[1],velocity[2])

from numpy.fft import rfftn, fftshift
fekin = fftshift(np.abs(rfftn(ekin)))**2

iws,nws,pws = shell_avg_3d(fekin)

iws = iws[nws != 0] / nws[nws != 0]
pws = pws[nws != 0] / nws[nws != 0]

for i,p in zip(iws,pws):
    print(i,p)
