#!/usr/bin/env python3

import numpy as np
from numpy.fft import rfftn, fftshift
import sys

import flash
from shellavg import shell_avg_3d
import ulz

try:
    fp = sys.argv[1]
except IndexError:
    print("usage: %s [checkpoint file] [nsamples (optional)]" % sys.argv[0])
    sys.exit(1)

try:
    nsamples = int(sys.argv[2])
except IndexError:
    nsamples = None

flash = flash.File(fp)

time = flash.realscalars['time']
step = flash.integerscalars['nstep']

c_s  = flash.realruntime['c_ambient']
rho0 = flash.realruntime['rho_ambient']

Vgrid   = np.prod(flash.gridsize)
Vcell   = np.prod(flash.cellsize) 
Vdomain = np.prod(flash.domainsize) 

density   = flash.data('dens')
velocity  = tuple(flash.data('vel'+dim) for dim in 'x y z'.split())

ekin  = 0.5*Vcell*density * ulz.norm(*velocity)
fekin = fftshift(np.abs(rfftn(ekin)))

radii,totals  = shell_avg_3d(fekin**2, nsamples)
powerspectrum = radii**2 * totals

np.savetxt(sys.stdout.buffer, np.array([radii,powerspectrum]).T)
