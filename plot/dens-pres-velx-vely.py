#!/usr/bin/env pyturbubox

import os
import sys
import numpy as np
import periodicbox, ulz # jmark
import pickle
import pathlib as pl
from collections import namedtuple
import multiprocessing as mpr
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 20})
from matplotlib import pyplot as plt

## ========================================================================= ##
## process commandline arguments

import argparse

pp = argparse.ArgumentParser(description = 'FLEXI Batch Plotter')

pp.add_argument(
    '--destdir',
    help='path to store: <dir>/%%03d.png',
    type=pl.Path, required=True,
)

pp.add_argument(
    '--title',
    type=str,
)

pp.add_argument(
    '--cachedir',
    help='path to cache min max calculations',
    type=pl.Path,
)

pp.add_argument(
    '--parallel',
    help='enable parallel processes: 0 --> max. n procs, > 0 --> set n procs',
    type=int,
    default=-1,
)

pp.add_argument(
    '--crosstime',
    help='crossing time scale: 1.0 (default)',
    type=float,
    default=1.0,
)

pp.add_argument(
    '--gather_min_max',
    help='flag to gather min max values beforehand',
    action='store_true',
)

pp.add_argument(
    '--ntasks',
    help='number of individual tasks',
    type=int,
)

pp.add_argument(
    '--crange',
    help='color range in the format "cdens = (-1,1), cmach = (-2,2)"',
)

pp.add_argument(
    'snapshots',
    help='list of snapshot files',
    type=pl.Path,nargs='*',
)

ARGV = pp.parse_args()

## ========================================================================= ##
## define tasks 

# interfaces used over subroutines below
Data   = namedtuple('Data', 'taskID time dyntime cmach cdens cpres cvort cekin')
CRange = namedtuple('CRange', 'cmach cdens cpres cvort cekin')

def calc_data(srcfp):
    box = periodicbox.File(srcfp.as_posix(), mode='r')

    dens, velx, vely, velz, pres = box.get_prims()

    fv = box.domainsize / np.array(dens.shape) # finite volume dimensions

    ax = 2
    cmach = np.nanmean(dens,axis=ax)
    cekin = np.nanmean(pres,axis=ax)
    cdens = np.nanmean(velx,axis=ax)
    cvort = np.nanmean(vely,axis=ax)
    cpres = cekin

    ax = 1
    cmach = np.nanmean(cmach,axis=ax)
    cekin = np.nanmean(cekin,axis=ax)
    cdens = np.nanmean(cdens,axis=ax)
    cvort = np.nanmean(cvort,axis=ax)
    cpres = cekin

    # ax = 2
    # cekin = np.log10(np.nanmean(ekin,axis=ax))
    # cmach = np.log10(np.nanmean(mach,axis=ax))
    # cdens = np.log10(np.nanmean(dens,axis=ax))
    # cpres = np.log10(np.nanmean(pres,axis=ax))
    # cvort = np.log10(np.nanmean(vort,axis=ax))

    # ax = ekin.shape[2] // 2
    # cekin = np.log10(ekin[:,:,ax])
    # cmach = np.log10(mach[:,:,ax])
    # cdens = np.log10(dens[:,:,ax])
    # cpres = np.log10(pres[:,:,ax])
    # cvort = np.log10(vort[:,:,ax])

    print('Finnished:  ', srcfp, flush=True)

    return Data(
        taskID   = -1,
        time     = box.time,
        dyntime  = box.time / ARGV.crosstime,
        cekin    = cekin,
        cmach    = cmach,
        cdens    = cdens,
        cpres    = cpres,
        cvort    = cvort
    )

def min_max(taskID, srcfp):
    data = calc_data(srcfp)

    result = Data(
        taskID = taskID, time = data.time, dyntime = data.dyntime,
        cekin  = ( np.min(data.cekin), np.max(data.cekin) ),
        cmach  = ( np.min(data.cmach), np.max(data.cmach) ),
        cdens  = ( np.min(data.cdens), np.max(data.cdens) ),
        cpres  = ( np.min(data.cpres), np.max(data.cpres) ),
        cvort  = ( np.min(data.cvort), np.max(data.cvort) ),
    )

    print(ulz.flatten_dict(result._asdict()), flush=True)
    
    return result

def mkplot(taskID, srcfp, crange):
    proc = mpr.current_process()
    data = calc_data(srcfp)

    subplt = [2,2,0]
    fig = plt.figure(figsize=(20, 18))

    frameid = taskID+1
    nframes = ARGV.ntasks if ARGV.ntasks else len(ARGV.snapshots)

    if ARGV.title is None:
        title = "dyntime: % 2.4f | frame: %03d/%03d" % (data.dyntime/ ARGV.crosstime, frameid, nframes)
    else:
        title = "%s (dyntime: % 2.4f | frame: %03d/%03d)" % (ARGV.title, data.dyntime/ ARGV.crosstime, frameid, nframes)

    plt.suptitle(title, fontsize='x-large').set_y(1.01)

    def plot(data, crange, title):
        subplt[2] += 1
        ax = fig.add_subplot(*subplt)
        ax.set_title(title)
        ax.set_xlabel('x index')
        ax.set_ylabel('y index')
        ax.grid()

        # img = ax.imshow(data,
        #     vmin = crange[0],
        #     vmax = crange[1],
        #     interpolation = 'none',
        #     cmap = plt.get_cmap('cubehelix'),
        # )

        #plt.colorbar(img,fraction=0.0456, pad=0.04, format='%1.2f')

        ax.plot(data, lw=2)

    plot(data.cmach, crange.cmach, title='column mach number (log10)')
    plot(data.cdens, crange.cdens, title='column density (log10)')
    #plot(data.cpres, crange.cpres, title='column pressure (log10)')
    plot(data.cekin, crange.cekin, title='column kin. energy (log10)')
    plot(data.cvort, crange.cvort, title='column vorticity (log10)')

    fig.tight_layout()

    plotfp = ARGV.destdir / srcfp.with_suffix('.png').name
    plt.savefig(str(plotfp), bbox_inches='tight')
    plt.close()
    print('Finnished:  ', str(plotfp), flush=True)

## ========================================================================= ##
## activate caching

if ARGV.cachedir:
    class mask: pass

    mask.calc_data = calc_data
    def calc_data(srcfp):
        cachefp = ARGV.cachedir / srcfp.with_suffix('.cdata.cache.pickle').name
        print('Processing: ', cachefp, flush=True)
        return ulz.cache(srcfp, cachefp, mask.calc_data, srcfp)

    mask.min_max = min_max
    def min_max(taskID, srcfp):
        cachefp = ARGV.cachedir / srcfp.with_suffix('.minmax.cache.pickle').name
        print('Processing: ', cachefp, flush=True)
        retval = ulz.cache(srcfp.as_posix(), cachefp.as_posix(), mask.min_max, taskID, srcfp)
        print('Finnished:  ', cachefp, flush=True)
        return retval

    mask.mkplot = mkplot
    def mkplot(taskID, srcfp, crange):
        plotfp = ARGV.destdir / srcfp.with_suffix('.png').name
        print('Processing: ', plotfp, flush=True)
        if plotfp.exists() and plotfp.stat().st_mtime > srcfp.stat().st_mtime:
            return
        return mask.mkplot(taskID, srcfp, crange)

## ========================================================================= ##
## set color range defaults

crange = {k: (None,None) for k in CRange._fields}

if ARGV.crange:
    # DANGER: using 'eval' on tainted data poses a security risk
    crange.update(eval('dict(%s)' % ARGV.crange))

crange = CRange(**crange)

## ========================================================================= ##
## gather minimun and maximum values

if ARGV.gather_min_max:
    if ARGV.parallel >= 0:
        def task(args):
            return min_max(*args)
        nprocs = None if ARGV.parallel == 0 else ARGV.parallel
        tmp = mp.Pool(nprocs).map(task,enumerate(ARGV.snapshots))
    else:
        tmp = [min_max(i,x) for i,x in enumerate(ARGV.snapshots)]

    def sanitize(dname, i):
        return [x for x in (getattr(X,dname)[i] for X in tmp) if not np.isnan(x) or np.isinf(x)]

    crange = CRange(
        cmach = ( np.min(sanitize('cmach',0)), np.max(sanitize('cmach',1)) ),
        cdens = ( np.min(sanitize('cdens',0)), np.max(sanitize('cdens',1)) ),
        cekin = ( np.min(sanitize('cekin',0)), np.max(sanitize('cekin',1)) ),
        cpres = ( np.min(sanitize('cpres',0)), np.max(sanitize('cpres',1)) ),
        cvort = ( np.min(sanitize('cvort',0)), np.max(sanitize('cvort',1)) ),
    )

## ========================================================================= ##
## do plotting

if ARGV.parallel >= 0:
    def task(args):
        return mkplot(args[0], args[1], crange)
    nprocs = None if ARGV.parallel == 0 else ARGV.parallel
    mpr.Pool(nprocs,maxtasksperchild=1).map(task,enumerate(ARGV.snapshots))
else:
    [mkplot(i,x,crange) for i,x in enumerate(ARGV.snapshots)]
