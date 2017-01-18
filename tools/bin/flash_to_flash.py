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
    zoomfactor = 1/2
    return scipy.ndimage.interpolation.zoom(src, zoomfactor, order=1, mode='wrap')

def scale(src):
    scalefactor = 1
    return scalefactor * src

def blur(src):
    blurfactor = 2
    return scipy.ndimage.filters.gaussian_filter(src, blurfactor)

# =========================================================================== #

with dslopts.Manager(scope=globals()) as mgr:
    mgr.add(name='srcfp', desc='source flash file path', type=dslopts.types.ExistingPath)
    mgr.add(name='snkfp', desc='  sink flash file path', type=dslopts.types.ExistingPath)

print("  var   |       min    ->    min     |       max    ->    max     ")
print("  ------|----------------------------|----------------------------")

with flash.File(srcfp,mode='r') as srcfls:
    with flash.File(snkfp,mode='r+') as snkfls:
        srcdens, srcvelx, srcvely, srcvelz, srcpres = srcfls.get_prims()
        srcener = srcfls.get_data('ener')

        gamma = srcfls.params['gamma']

        snkdens = zoom(srcdens)
        snkvelx = zoom(srcvelx)
        snkvely = zoom(srcvely)
        snkvelz = zoom(srcvelz)
        snkpres = zoom(srcpres)
        snkener = snkpres/(gamma-1)/snkdens + snkdens/2 * (snkvelx**2 + snkvely**2 + snkvelz**2)

        srcs = [srcdens, srcvelx, srcvely, srcvelz, srcpres, srcener]
        snks = [snkdens, snkvelx, snkvely, snkvelz, snkpres, snkener]

        for dbname, src, snk in zip('dens velx vely velz pres ener'.split(), srcs, snks):
            print("  %s  | % 12.5f % 12.5f  | % 12.5f % 12.5f" % (dbname, src.min(), snk.min(), src.max(), snk.max()))
            snkfls.set_data(dbname, snk)
