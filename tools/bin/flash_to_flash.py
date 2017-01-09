#!/usr/bin/env python3

import time
import flash
import hopr
import flexi
import gausslobatto
import scipy.interpolate
from scipy import ndimage
import numpy as np
import sys
import ulz
import interpolate
import dslopts
import pathlib
import scipy.ndimage
               
# =========================================================================== #

def ExistingPath(arg):
    pth = pathlib.Path(arg)
    if not pth.exists():
        raise OSError("'%s' does not exists!" % pth)
    return arg

with dslopts.Manager(scope=globals()) as mgr:
    mgr.add(name='srcfp', desc='source flash file path', type=ExistingPath)
    mgr.add(name='snkfp', desc='  sink flash file path', type=ExistingPath)

print("  var   |       min    ->    min     |       max    ->    max     ")
print("  ------|----------------------------|----------------------------")

blurfactor = 2
scalefactor = 0.2

with flash.File(srcfp,mode='r') as srcfls:
    with flash.File(snkfp,mode='r+') as snkfls:
        for dbname in 'dens velx vely velz pres'.split():
            src = srcfls.get_data(dbname)
            src = scipy.ndimage.interpolation.zoom(src, 0.5,order=1, mode='wrap')

            if dbname in "dens pres":
                snk = np.ones_like(src)
                
            if dbname in "velx vely velz":
                snk = scalefactor * scipy.ndimage.filters.gaussian_filter(src, blurfactor)

            print("  %s  | % 12.5f % 12.5f  | % 12.5f % 12.5f" % (dbname, src.min(), snk.min(), src.max(), snk.max()))
            
            snkfls.set_data(dbname, snk)
