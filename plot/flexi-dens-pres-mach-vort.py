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
    'snapshots',
    help='list of snapshot files',
    type=str,nargs='*',
)

pp.add_argument('--parallel', help='parallel processing', action='store_true')

cmdargs = pp.parse_args()

## ========================================================================= ##
## define tasks 

def calc_data(srcfp):

    box = flexi.PeriodicBox(srcfp)

    time = box.time
    c_s  = 1.
    MACH = 8.
    turntime = time / (np.mean(box.domainsize) / c_s / MACH)

    #dens, velx, vely, velz, pres = box.get_prims()
    dens, velx, vely, velz, pres = box.get_prims_fv()

    dens[dens < 0] = 1e-5
    pres[pres < 0] = 1e-5

    fv = box.domainsize / np.array(dens.shape) # finite volume dimensions

    mach = np.sqrt(velx**2+vely**2+velz**2)/np.sqrt(pres/dens)
    vort = np.mean(fv)**5/12.0 * np.abs(dens) * ulz.norm(*ulz.curl(velx,vely,velz,fv[0],fv[1],fv[2]))

    mach[np.isinf(mach)] = np.nan
    dens[np.isinf(dens)] = np.nan
    pres[np.isinf(pres)] = np.nan
    vort[np.isinf(vort)] = np.nan

    ax = 2
    cmach = np.log10(np.nanmean(mach,axis=ax))
    cdens = np.log10(np.nanmean(dens,axis=ax))
    cpres = np.log10(np.nanmean(pres,axis=ax))
    cvort = np.log10(np.nanmean(vort,axis=ax))

    return {
        'time':     time,
        'turntime': turntime,
        'cmach':    cmach,
        'cdens':    cdens,
        'cpres':    cpres,
        'cvort':    cvort
    }

def min_max(taskID, srcfp):

    data = calc_data(srcfp)

    result = [ 
        taskID,data['time'], 

        np.min(data['cmach']),
        np.max(data['cmach']),
 
        np.min(data['cdens']),
        np.max(data['cdens']),

        np.min(data['cpres']),
        np.max(data['cpres']),

        np.min(data['cvort']),
        np.max(data['cvort']),
    ]

    print(*result, flush=True)
    
    return result

def mkplot(taskID, srcfp):
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams.update({'font.size': 20})
    from matplotlib import pyplot as plt

    data = calc_data(srcfp)

    time = data['time']
    turntime = data['turntime']

    subplt = [2,2,0]
    fig = plt.figure(figsize=(20, 18))

    if cmdargs.title is None:
        title = "dyntime: % 2.4f | frame: %03d/%03d" % (cmdargs.title, turntime, taskID+1, len(cmdargs.snapshots))
    else:
        title = "%s (dyntime: % 2.4f | frame: %03d/%03d)" % (cmdargs.title, turntime, taskID+1, len(cmdargs.snapshots))

    plt.suptitle(title, fontsize='x-large').set_y(1.01)

    def plot(data, title, crange=None):
        subplt[2] += 1
        ax = fig.add_subplot(*subplt)
        ax.set_title(title)
        ax.set_xlabel('x index'); ax.set_ylabel('y index')

        p = {'cmap': plt.get_cmap('cubehelix'), 'interpolation': 'none'}
        if crange is None:
            img = ax.imshow(data, **p)
        else:
            img = ax.imshow(data, vmin=crange[0], vmax=crange[1], **p)
        plt.colorbar(img,fraction=0.0456, pad=0.04, format='%1.2f')

    plot(data['cmach'], 'column mach number (log10)', crange['cmach'])
    plot(data['cdens'], 'column density (log10)', crange['cdens'])
    plot(data['cpres'], 'column pressure (log10)', crange['cpres'])
    plot(data['cvort'], 'column vorticity (log10)', crange['cvort'])

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

def task(args):
    return min_max(*args)

if cmdargs.parallel:
    import multiprocessing as mp
    tmp = mp.Pool().map(task,enumerate(cmdargs.snapshots))
else:
    tmp = [task(x) for x in enumerate(cmdargs.snapshots)]

tmp = np.array(tmp)
i = ulz.mkincr(start=2)

crange = {
    'cmach': (np.min(tmp[:,next(i)]), np.max(tmp[:,next(i)])),
    'cdens': (np.min(tmp[:,next(i)]), np.max(tmp[:,next(i)])),
    'cpres': (np.min(tmp[:,next(i)]), np.max(tmp[:,next(i)])),
    'cvort': (np.min(tmp[:,next(i)]), np.max(tmp[:,next(i)]))
}

## ========================================================================= ##
## do plotting

def task(args):
    return mkplot(*args)

if cmdargs.parallel:
    import multiprocessing as mpr
    mpr.Pool().map(task,enumerate(cmdargs.snapshots))
else:
    [task(x) for x in enumerate(cmdargs.snapshots)]
