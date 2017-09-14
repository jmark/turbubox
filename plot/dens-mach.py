#!/usr/bin/env pyturbubox

import os
import box as dbox
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
Data   = namedtuple('Data', 'taskID time dyntime cdens cmach')
CRange = namedtuple('CRange', 'cdens cico cipc cmach')

def calc_data(srcfp):
    box = periodicbox.File(srcfp.as_posix(), mode='r')

    dens, velx, vely, velz, pres = box.get_prims()

    ico = box.get_data('ih2 ')
    icp = box.get_data('tdus')

    dens[dens < 0] = 1e-5
    pres[pres < 0] = 1e-5

    mach = np.sqrt(velx**2+vely**2+velz**2)/np.sqrt(pres/dens)

    ax = 2
    return dbox.Box(
        taskID   = -1,
        time     = box.time,
        dyntime  = box.time / ARGV.crosstime,
        cdens    = np.log10(np.sum(dens,axis=ax)),
        cmach    = np.log10(np.mean(mach,axis=ax)),
        #cico     = np.log10(np.sum(ico,axis=ax)),
        #cicp     = np.log10(np.sum(ico,axis=ax)),
        cico     = ico[:,:,64],
        cicp     = icp[:,:,64],
    )

def min_max(taskID, srcfp):
    data = calc_data(srcfp)

    result = Data(
        taskID = taskID, time = data.time, dyntime = data.dyntime,
        cdens  = ( np.min(data.cdens), np.max(data.cdens) ),
        cmach  = ( np.min(data.cmach), np.max(data.cmach) ),
    )

    print(ulz.flatten_dict(result._asdict()), flush=True)
    
    return result

def mkplot(taskID, srcfp, crange):
    proc = mpr.current_process()
    data = calc_data(srcfp)
    print('Finnished:  ', srcfp, flush=True)

    subplt = [1,2,0]
    fig = plt.figure(figsize=(20, 9.5))

    frameid = taskID+1
    nframes = ARGV.ntasks if ARGV.ntasks else len(ARGV.snapshots)

    if ARGV.title is None:
        title = "simtime: {:1.4e} | frame: {:03d}/{:03d}".format(data.time, frameid, nframes)
        #title = "dyntime: % 2.4f | frame: %03d/%03d" % (data.dyntime/ ARGV.crosstime, frameid, nframes)
    else:
        title = "%s (dyntime: % 2.4f | frame: %03d/%03d)" % (ARGV.title, data.dyntime, frameid, nframes)

    plt.suptitle(title, fontsize='x-large').set_y(1.01)

    def plot(data, crange, title):
        subplt[2] += 1
        ax = fig.add_subplot(*subplt)
        ax.set_title(title)
        ax.set_xlabel('x index')
        ax.set_ylabel('y index')

        img = ax.imshow(data,
            #vmin = crange[0],
            #vmax = crange[1],
            interpolation = 'none',
            cmap = plt.get_cmap('cubehelix'),
        )

        plt.colorbar(img,fraction=0.0456, pad=0.04, format='%1.2f')

    #plot(data.cdens, crange.cdens, title='column sum density (log10)')
    plot(data.cico, crange.cico, title='column sum ico (log10)')
    plot(data.cico, crange.cico, title='column sum icp (log10)')
    #plot(data.cmach, crange.cmach, title='column mean mach number (log10)')

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
        cdens = ( np.min(sanitize('cdens',0)), np.max(sanitize('cdens',1)) ),
        cmach = ( np.min(sanitize('cmach',0)), np.max(sanitize('cmach',1)) ),
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
