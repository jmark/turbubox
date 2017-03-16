#!/usr/bin/env pyturbubox

import os
import sys
import numpy as np
import flexi, ulz # jmark
import pickle

## ========================================================================= ##
## process commandline arguments

import argparse

pp = argparse.ArgumentParser(description = 'FLEXI Batch Plotter')

pp.add_argument(
    '--destdir',
    help='path to store: <dir>/%%03d.png',
    type=str, required=True,
)

pp.add_argument('--title',type=str,)

pp.add_argument(
    '--cachedir',
    help='path to cache min max calculations',
    type=str,
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
    'snapshots',
    help='list of snapshot files',
    type=str,nargs='*',
)

cmdargs = pp.parse_args()

## ========================================================================= ##
## define tasks 

def calc_data(srcfp):
    box = flexi.PeriodicBox(srcfp)

    #dens, velx, vely, velz, pres = box.get_prims()
    dens, velx, vely, velz, pres = box.get_prims_fv()

    dens[dens < 0] = 1e-5
    pres[pres < 0] = 1e-5

    fv = box.domainsize / np.array(dens.shape) # finite volume dimensions

    mach = np.sqrt(velx**2+vely**2+velz**2)/np.sqrt(pres/dens)
    vort = np.mean(fv)**5/12.0 * np.abs(dens) * ulz.norm(*ulz.curl(velx,vely,velz,fv[0],fv[1],fv[2]))

    ax = 2
    cmach = np.log10(np.nanmean(mach,axis=ax))
    cdens = np.log10(np.nanmean(dens,axis=ax))
    cpres = np.log10(np.nanmean(pres,axis=ax))
    cvort = np.log10(np.nanmean(vort,axis=ax))

    return {
        'time':     box.time,
        'dyntime':  box.time / cmdargs.crosstime,
        'cmach':    cmach,
        'cdens':    cdens,
        'cpres':    cpres,
        'cvort':    cvort
    }

def min_max(taskID, srcfp):
    data = calc_data(srcfp)

    result = {
        'taskID': taskID, 'time': data['time'], 
        'cmach': (np.min(data['cmach']), np.max(data['cmach'])),
        'cdens': (np.min(data['cdens']), np.max(data['cdens'])),
        'cpres': (np.min(data['cpres']), np.max(data['cpres'])),
        'cvort': (np.min(data['cvort']), np.max(data['cvort'])),
    }

    print(ulz.flatten_dict(result), flush=True)
    
    return result

def mkplot(taskID, srcfp):
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams.update({'font.size': 20})
    from matplotlib import pyplot as plt

    data = calc_data(srcfp)

    time = data['time']
    dyntime = data['dyntime']

    subplt = [2,2,0]
    fig = plt.figure(figsize=(20, 18))

    if cmdargs.title is None:
        title = "dyntime: % 2.4f | frame: %03d/%03d" % (dyntime, taskID+1, len(cmdargs.snapshots))
    else:
        title = "%s (dyntime: % 2.4f | frame: %03d/%03d)" % (cmdargs.title, dyntime, taskID+1, len(cmdargs.snapshots))

    plt.suptitle(title, fontsize='x-large').set_y(1.01)

    def plot(dname, title):
        subplt[2] += 1
        ax = fig.add_subplot(*subplt)
        ax.set_title(title)
        ax.set_xlabel('x index')
        ax.set_ylabel('y index')

        p = {
            'cmap': plt.get_cmap('cubehelix'),
            'vmin': crange[dname][0] if crange is not None else None, # 'crange' is defined further below
            'vmax': crange[dname][1] if crange is not None else None,
            'interpolation': 'none',
        }

        img = ax.imshow(data[dname], **p)
        plt.colorbar(img,fraction=0.0456, pad=0.04, format='%1.2f')

    plot('cmach', 'column mach number (log10)')
    plot('cdens', 'column density (log10)')
    plot('cpres', 'column pressure (log10)')
    plot('cvort', 'column vorticity (log10)')

    fig.tight_layout()

    plotfp = cmdargs.destdir + '/' + srcfp + '.png' 
    plt.savefig(plotfp,bbox_inches='tight')
    print(plotfp, flush=True)

    plt.close()

## ========================================================================= ##
## activate caching

if cmdargs.cachedir:
    class mask: pass

    mask.calc_data = calc_data
    def calc_data(srcfp):
        cachefp = cmdargs.cachedir + '/' + srcfp + '.cdata.cache.pickle'
        return ulz.cache(srcfp, cachefp, mask.calc_data, srcfp)

    mask.min_max = min_max
    def min_max(taskID, srcfp):
        cachefp = cmdargs.cachedir + '/' + srcfp + '.minmax.cache.pickle'
        return ulz.cache(srcfp, cachefp, mask.min_max, taskID, srcfp)

    mask.mkplot = mkplot
    def mkplot(taskID, srcfp):
        plotfp = cmdargs.destdir + '/' + srcfp + '.png' 
        if os.path.exists(plotfp) and os.path.getmtime(plotfp) > os.path.getmtime(srcfp):
            return
        return mask.mkplot(taskID, srcfp)

## ========================================================================= ##
## gather minimun and maximum values

crange = None

if cmdargs.gather_min_max:
    def task(args):
        return min_max(*args)

    if cmdargs.parallel >= 0:
        import multiprocessing as mp
        nprocs = None if cmdargs.parallel == 0 else cmdargs.parallel
        tmp = mp.Pool(nprocs).map(task,enumerate(cmdargs.snapshots))
    else:
        tmp = [task(x) for x in enumerate(cmdargs.snapshots)]

    def sanitize(dname, i):
        return [x for x in (X[dname][i] for X in tmp) if not (np.isnan(x) or np.isinf(x))]

    crange = {
        'cmach': ( np.min(sanitize('cmach',0)), np.max(sanitize('cmach',1)) ),
        'cdens': ( np.min(sanitize('cdens',0)), np.max(sanitize('cdens',1)) ),
        'cpres': ( np.min(sanitize('cpres',0)), np.max(sanitize('cpres',1)) ),
        'cvort': ( np.min(sanitize('cvort',0)), np.max(sanitize('cvort',1)) ),
    }

## ========================================================================= ##
## do plotting

def task(args):
    return mkplot(*args)

if cmdargs.parallel >= 0:
    import multiprocessing as mpr
    nprocs = None if cmdargs.parallel == 0 else cmdargs.parallel
    mpr.Pool(nprocs).map(task,enumerate(cmdargs.snapshots))
else:
    [task(x) for x in enumerate(cmdargs.snapshots)]
