#!/usr/bin/env python3

import time
import flash
import numpy as np
import sys
import ulz
import dslopts

import scipy.interpolate
import scipy.ndimage

# =========================================================================== #

def zoom(src):
    zoom = 2
    return scipy.ndimage.interpolation.zoom(src, zoom, order=1, mode='wrap')

def scale(src):
    scalefactor = 1.0
    return scalefactor * src

def blur(src)
    blurfactor = 2
    return scipy.ndimage.filters.gaussian_filter(src, blurfactor)

# =========================================================================== #

with dslopts.Manager(scope=globals()) as mgr:
    mgr.add(name='srcfp', desc='source flash file path', type=dslopts.types.ExistingPath)
    mgr.add(name='snkfp', desc='  sink flash file path', type=dslopts.types.ExistingPath)

print("  var   |       min    ->    min     |       max    ->    max     ")
print("  ------|----------------------------|----------------------------")

zoom = 2
blurfactor = 2
scalefactor = 1.0

with flash.File(srcfp,mode='r') as srcfls:
    with flash.File(snkfp,mode='r+') as snkfls:
        for dbname in 'dens velx vely velz pres'.split():
            src = srcfls.get_data(dbname)

            tmp = src
            tmp = zoom(tmp)

            # if dbname in "dens pres":
            #     tmp = np.ones_like(tmp)
            #     
            # if dbname in "velx vely velz":
            #     tmp = scale(tmp)

            snk = tmp
            print("  %s  | % 12.5f % 12.5f  | % 12.5f % 12.5f" % (dbname, src.min(), snk.min(), src.max(), snk.max()))
            snkfls.set_data(dbname, snk)
