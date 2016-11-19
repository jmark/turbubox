#!/usr/bin/env python3

import numpy as np
from numpy.fft import rfftn, fftshift
import sys
import pathlib

import flash
from shellavg import shell_avg_3d
import ulz
import dslopts

def ExistingPath(arg):
    pth = pathlib.Path(arg)
    if not pth.exists():
        raise OSError("'%s' does not exists!" % pth)
    return pth

def PositiveInt(arg):
    x = int(arg)
    if x >= 0:
        return x
    raise ValueError("'%d' must be positive!" % x)

with dslopts.Manager(scope=locals()) as mgr:
    mgr.add('fp'       ,'flash file path' ,ExistingPath)
    mgr.add('nsamples' ,'no. of samples'  ,PositiveInt, 0)

flash = flash.File(str(fp))

time = flash.realscalars['time']
step = flash.integerscalars['nstep']

c_s  = flash.realruntime['c_ambient']
rho0 = flash.realruntime['rho_ambient']

# c_s  = flash.realruntime['sim_cambient']
# rho0 = flash.realruntime['sim_rhoambient']

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
