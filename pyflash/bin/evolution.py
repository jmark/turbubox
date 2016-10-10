#!py.sh

import numpy as np
import sys

import flash
import vectoranalysis3D as va3D

def product(it):
    p = 1
    for i in it:
        p *= i
    return p

#print("# ","\t".join("count nstep time ekin evortex".split()))

for count,filepath in enumerate(sys.stdin):
    flash = flash.File(filepath.rstrip())

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
    vorticity = va3D.curl(velocity[0],velocity[1],velocity[2],CS[0],CS[1],CS[2])

    ekintotal = Vcell/2.0     * np.sum(density * va3D.norm(velocity[0],velocity[1],velocity[2]))
    evortotal = CS[0]**5/12.0 * np.sum(density * va3D.norm(vorticity[0],vorticity[1],vorticity[2]))

    print("\t".join(map(str,[count,step,time,ekintotal,evortotal])))
