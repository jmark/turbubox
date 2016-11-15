#!/usr/bin/env python3

from matplotlib import pylab as plt
from mpl_toolkits.mplot3d import axes3d
import sys

import numpy as np
from numpy.fft import rfftn, fftshift

import flash as FLASH
from shellavg import shell_avg_3d
import ulz


sys.argv.reverse()
progpath = sys.argv.pop()
flsfp = sys.argv.pop()
flsfp2 = sys.argv.pop()
sinkfp = sys.argv.pop()

flash = FLASH.File(flsfp)

time = flash.realscalars['time']
#step = flash.integerscalars['nstep']

c_s  = flash.realruntime['c_ambient']
rho0 = flash.realruntime['rho_ambient']

# c_s  = flash.realruntime['sim_cambient']
# rho0 = flash.realruntime['sim_rhoambient']

Vgrid   = np.prod(flash.gridsize)
Vcell   = np.prod(flash.cellsize) 
Vdomain = np.prod(flash.domainsize) 

density   = flash.data('dens')
velocity  = tuple(flash.data('vel'+dim) for dim in 'x y z'.split())

#ekin  = 0.5*Vcell*density * ulz.norm(*velocity)
#fekin = fftshift(np.abs(rfftn(ekin)))

#input  = density**(1/3) * np.sqrt(ulz.norm(*velocity))
input  = np.sqrt(ulz.norm(*velocity))
#vels  = ulz.norm(*velocity)
fvels = fftshift(np.abs(rfftn(input)))

nsamples = 200
radii,totals = shell_avg_3d(fvels**2, nsamples)

rs = radii[1:]
ps = rs**2 * totals[1:]


plt.grid()
plt.xlim(1,1000)

xs = rs
ys = ps/rs**(-5/3)
plt.loglog(xs,ys, '-', label='b3')


#########################

flash = FLASH.File(flsfp2)

time = flash.realscalars['time']
#step = flash.integerscalars['nstep']

c_s  = flash.realruntime['c_ambient']
rho0 = flash.realruntime['rho_ambient']

# c_s  = flash.realruntime['sim_cambient']
# rho0 = flash.realruntime['sim_rhoambient']

Vgrid   = np.prod(flash.gridsize)
Vcell   = np.prod(flash.cellsize) 
Vdomain = np.prod(flash.domainsize) 

density   = flash.data('dens')
velocity  = tuple(flash.data('vel'+dim) for dim in 'x y z'.split())

#ekin  = 0.5*Vcell*density * ulz.norm(*velocity)
#fekin = fftshift(np.abs(rfftn(ekin)))

#input  = density**(1/3) * np.sqrt(ulz.norm(*velocity))
input  = np.sqrt(ulz.norm(*velocity))
#vels  = ulz.norm(*velocity)
fvels = fftshift(np.abs(rfftn(input)))

nsamples = 200
radii,totals = shell_avg_3d(fvels**2, nsamples)

rs = radii[1:]
ps = rs**2 * totals[1:]

#########################

xs = rs
ys = ps/rs**(-2)
plt.loglog(xs,ys, '-', label='b5')

plt.legend(loc='upper right')
plt.savefig(sinkfp,bbox_inches='tight')
