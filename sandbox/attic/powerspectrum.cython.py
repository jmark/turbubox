#!/usr/bin/env python3

import numpy as np
import sys

import flash
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

flash = flash.File(fp)

time = flash.realscalars['time']
step = flash.integerscalars['nstep']

c_s  = flash.realruntime['c_ambient']
rho0 = flash.realruntime['rho_ambient']

GS = flash.gridsize
CS = flash.cellsize
DS = flash.domainsize

Vgrid   = product(GS)
Vcell   = product(CS) 
Vdomain = product(DS) 

density   = flash.data('dens')
velocity  = [flash.data('vel'+dim) for dim in 'x y z'.split()]

ekin = 0.5*Vcell*density * va3D.norm(velocity[0],velocity[1],velocity[2])

from numpy.fft import rfftn, fftshift
fekin = fftshift(np.abs(rfftn(ekin)))**2

iws,nws,pws = shell_avg_3d(fekin)

iws = iws[nws != 0] / nws[nws != 0]
pws = pws[nws != 0] / nws[nws != 0]

for i,p in zip(iws,pws):
    print(i,p)
