#!/usr/bin/env python3

# stdlib
import sys
import numpy as np
from matplotlib import pyplot as plt

# jmark
import flash, ulz, dslopts

with dslopts.Manager(scope=globals()) as mgr:
    pass
    #mgr.add('flashfilepath')
    #mgr.add('flashfilepath')

def plot(flashfilepath):
    fls = flash.File(flashfilepath)

    mach = 10
    time = fls.realscalars['time']
    step = fls.integerscalars['nstep']
    snds = fls.realruntime['c_ambient']
    dns0 = fls.realruntime['rho_ambient']
    dyns = snds * mach # dynamic speed
    dynt = fls.domainsize[0] / dyns # dynamic time scale
    turn = time / dynt # turn time
    cell = np.prod(fls.cellsize) # cell volume

    dens = fls.data('dens')
    pres = fls.data('pres')
    velx = fls.data('velx')
    vely = fls.data('vely')
    velz = fls.data('velz')

    ekin = cell/2 * dens * (velx**2+vely**2+velz**2)

    data = pres
    dist, bins = np.histogram(data,
        bins=10000,
        #range=(0,4),
        density=True)

    xs = bins[0:-1] + np.abs(bins[1]-bins[0])/2
    ys = dist

    xs = np.log10(xs)
    ys = np.log10(ys)

    plt.plot(xs,ys)

fps = _ignored_

for fp in fps:
    plot(fp)

plt.show()
