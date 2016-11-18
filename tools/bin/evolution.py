#!/usr/bin/env python3

import numpy as np
import sys

import flash as fl
import ulz

for count,filepath in enumerate(sys.stdin):
    flash = fl.File(filepath.rstrip())

    time = flash.realscalars['time']
    step = flash.integerscalars['nstep']

    c_s  = flash.realruntime['c_ambient']
    rho0 = flash.realruntime['rho_ambient']

    GS = flash.gridsize
    DS = flash.domainsize
    CS = flash.cellsize

    Vgrid   = np.prod(GS)
    Vcell   = np.prod(CS) 
    Vdomain = np.prod(DS) 

    density   = flash.data('dens')
    velocity  = [flash.data('vel'+dim) for dim in 'x y z'.split()]
    mach      = np.sqrt(np.sum(density * ulz.norm(*velocity))/np.sum(density)) / c_s

    ekintotal = Vcell/2.0 * np.sum(density * ulz.norm(*velocity))

    vorticity = ulz.curl(*velocity,*tuple(CS))
    evortotal = CS[0]**5/12.0 * np.sum(density * ulz.norm(*vorticity))

    print("\t".join(map(str,[count,step,time,mach,ekintotal,evortotal])), flush=True)
