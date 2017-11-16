#!/usr/bin/env pyturbubox

import os
import sys
import numpy as np
import cubicle, ulz # jmark
import pickle
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 10})
from matplotlib import pyplot as plt

import box as bbox

## ========================================================================= ##
## process commandline arguments

import argparse

pp = argparse.ArgumentParser(description = 'FLEXI Batch Plotter')

pp.add_argument(
    '--destdir',
    help='path to store: <dir>/%%03d.png',
    type=Path, required=True,
)

pp.add_argument(
    '--profile',
    type=str,
    required=True,
)

pp.add_argument(
    '--skip',
    type=bool,
    default=False,
    help='skip already produced files',
)

pp.add_argument(
    '--parallel',
    help='enable parallel processes: 0 --> max. n procs, > 0 --> set n procs',
    type=int,
    default=-1,
)

pp.add_argument(
    '--snapshots',
    help='list of snapshot files',
    type=Path, nargs='*', required=True
)

ARGV = pp.parse_args()

## ========================================================================= ##
## define tasks 

if ARGV.profile == 'density':
    def calc_data(srcfp):
        box = cubicle.File(srcfp.as_posix(), mode='r')

        dens = box.get_data('dens')

        return bbox.Box(
            time     = box.time,
            data     = np.log10(np.nansum(dens,axis=2)),
            #crange   = (-26,25.5),
            crange   = (None,None),
            title    = 'simtime: {time:1.2e} | frame: {frameid}/{nframes}',
            subtitle = 'column total density (log10)',
        )

elif ARGV.profile == 'temperature-slice':
    def calc_data(srcfp):
        box = cubicle.File(srcfp.as_posix(), mode='r')

        dens, velx, vely, velz, pres = box.get_prims()

        #for k,v in box.realruntime.items():
        #    print(k, ' => ', v)

        temp = pres/dens

        return bbox.Box(
            time     = box.time,
            data     = np.log10(temp[:,:,0]),
            #data     = np.log10(temp[:,:,temp.shape[2]//2]),
            #crange   = (-26,25.5),
            crange   = (None,None),
            title    = 'simtime: {time:1.2e} | frame: {frameid}/{nframes}',
            subtitle = 'temperature (R*pres/dens) slice at z = 0 (log10)',
        )

elif ARGV.profile == 'mach':
    def calc_data(srcfp):
        box = cubicle.File(srcfp.as_posix(), mode='r')

        dens, velx, vely, velz, pres = box.get_prims()

        dens[dens < 0] = 1e-5
        pres[pres < 0] = 1e-5

        mach = np.sqrt(velx**2+vely**2+velz**2)/np.sqrt(pres/dens)

        return bbox.Box(
            time     = box.time,
            data     = np.log10(np.mean(mach,axis=2)),
            #crange   = (0.0,1.0),
            crange   = (None,None),
            title    = 'simtime: {time:1.2e} | frame: {frameid}/{nframes}',
            subtitle = 'column mean Mach (log10)',
        )

elif ARGV.profile in 'ico icp ih2 iha ihp':
    def calc_data(srcfp):
        box = cubicle.File(srcfp.as_posix(), mode='r')

        data = box.get_data(ARGV.profile + ' ')

        return bbox.Box(
            time     = box.time,
            data     = np.log10(np.sum(data,axis=2)),
            #crange   = (0.0,1.0),
            crange   = (None,None),
            title    = 'simtime: {time:1.2e} | frame: {frameid}/{nframes}',
            subtitle = 'column total {} (log10)'.format(ARGV.profile.upper()),
        )

else:
    raise NotImplementedError('Unknown profile: {}'.format(ARGV.profile))

def task(args):
    taskID, srcfp = args

    snkfp = ARGV.destdir / srcfp.with_suffix('.png').name
    if ARGV.skip and snkfp.exists() and snkfp.stat().st_mtime > srcfp.stat().st_mtime: return

    data = calc_data(srcfp)

    subplt = [1,1,1]
    fig = plt.figure(figsize=(10, 9.5))

    data.frameid = taskID+1
    data.nframes = len(ARGV.snapshots)

    plt.suptitle(data.title.format(**data), fontsize='large').set_y(1.00)

    ax = fig.add_subplot(*subplt)
    ax.set_title(data.subtitle)
    ax.set_xlabel('x index')
    ax.set_ylabel('y index')

    img = ax.imshow(data.data,
        vmin = data.crange[0],
        vmax = data.crange[1],
        interpolation = 'none',
        cmap = plt.get_cmap('cubehelix'),
    )

    plt.colorbar(img,fraction=0.0456, pad=0.04, format='%1.2f')

    fig.tight_layout()

    plt.savefig(str(snkfp), bbox_inches='tight')
    plt.close()
    print('Finnished:  ', str(snkfp), flush=True)

# ========================================================================= ##
## do plotting

if ARGV.parallel < 0:
    for i,x in enumerate(ARGV.snapshots): task((i,x))
else:
    import multiprocessing as mpr
    nprocs = None if ARGV.parallel == 0 else ARGV.parallel
    mpr.Pool(nprocs,maxtasksperchild=1).map(task,enumerate(ARGV.snapshots))
