#!py.sh

import numpy as np
import sys

from pyflash import FlashFile
import vectoranalysis3D as va3D

def product(it):
    p = 1
    for i in it:
        p *= i
    return p

print("# ","\t".join("count nstep time ekin evortex".split()))

for count,filepath in enumerate(sys.stdin):

    flash = FlashFile(filepath.rstrip())

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
    vorticity = va3D.curl(velocity[0],velocity[1],velocity[2],CS[0],CS[1],CS[2])

    ekintotal = Vcell/2.0     * np.sum(density * va3D.norm(velocity[0],velocity[1],velocity[2]))
    evortotal = CS[0]**5/12.0 * np.sum(density * va3D.norm(vorticity[0],vorticity[1],vorticity[2]))

    out = [count,step,time,ekintotal,evortotal]
    print("  ","\t".join(map(str,out)))
