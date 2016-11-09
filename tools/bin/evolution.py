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
    CS = flash.cellsize
    DS = flash.domainsize

    Vgrid   = np.prod(GS)
    Vcell   = np.prod(CS) 
    Vdomain = np.prod(DS) 

    density   = flash.data('dens')
    velocity  = [flash.data('vel'+dim) for dim in 'x y z'.split()]
    vorticity = ulz.curl(velocity[0],velocity[1],velocity[2],CS[0],CS[1],CS[2])

    ekintotal = Vcell/2.0     * np.sum(density * ulz.norm(velocity[0],velocity[1],velocity[2]))
    evortotal = CS[0]**5/12.0 * np.sum(density * ulz.norm(vorticity[0],vorticity[1],vorticity[2]))

    print("\t".join(map(str,[count,step,time,ekintotal,evortotal])), flush=True)