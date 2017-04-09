#!/usr/bin/env pyturbubox

import time
import flash
import numpy as np
import sys
import ulz
import dslopts

import scipy.interpolate
import scipy.ndimage

# =========================================================================== #

def pipe(source, *filters):
    sink = source
    for filter in filters:
        sink = filter(sink) 
    return sink

def blur(src):
    factor = 2
    return scipy.ndimage.filters.gaussian_filter(src, factor)

def zoom(src):
    factor = 1
    #factor = 1/2
    #factor = 396 / 256
    return scipy.ndimage.interpolation.zoom(src, factor, order=1, mode='wrap')

def scale(src):
    factor = 5
    return factor * src

def transform_normal(box):
    avg = box.mean()
    return avg / np.mean(box) * pipe(box,blur,zoom) 

def transform_scaled(box):
    return pipe(box,blur,scale,zoom)

# =========================================================================== #

with dslopts.Manager(scope=globals()) as mgr:
    mgr.add(name='srcfp', desc='source flash file path', type=dslopts.types.ExistingPath)
    mgr.add(name='snkfp', desc='  sink flash file path', type=dslopts.types.ExistingPath)

print("  var   |       min    ->    min     |       max    ->    max     ")
print("  ------|----------------------------|----------------------------")

with flash.File(srcfp,mode='r') as srcfls:
    with flash.File(snkfp,mode='r+') as snkfls:
        srcdens, srcvelx, srcvely, srcvelz, srcpres = srcfls.get_prims()
        srceint = srcfls.get_data('eint')
        srcener = srcfls.get_data('ener')

        gamma = srcfls.params['gamma']

        snkdens = transform_normal(srcdens)
        snkvelx = transform_scaled(srcvelx)
        snkvely = transform_scaled(srcvely)
        snkvelz = transform_scaled(srcvelz)
        snkpres = transform_normal(srcpres)

        snkeint = snkpres/(gamma-1)/snkdens
        snkener = snkpres/(gamma-1)/snkdens + 0.5*(snkvelx**2 + snkvely**2 + snkvelz**2)

        srcs = [srcdens, srcvelx, srcvely, srcvelz, srcpres, srceint, srcener]
        snks = [snkdens, snkvelx, snkvely, snkvelz, snkpres, snkeint, snkener]

        for dbname, src, snk in zip('dens velx vely velz pres eint ener'.split(), srcs, snks):
            print("  %s  | % 12.5f % 12.5f  | % 12.5f % 12.5f" % (dbname, src.min(), snk.min(), src.max(), snk.max()))
            snkfls.set_data(dbname, snk)
