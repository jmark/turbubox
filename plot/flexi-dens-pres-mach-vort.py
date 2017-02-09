#!/usr/bin/env python3

# stdlib
import os
import sys
import numpy as np
import multiprocessing as mpr

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

matplotlib.rcParams.update({'font.size': 20})

# jmark
import flexi, ulz, dslopts

SOLVER = 'DG'
MACH = 2

def calc_data(srcfp):
    box = flexi.PeriodicBox(srcfp)

    time = box.time
    c_s  = 1
    turntime = time / (np.mean(box.domainsize) / c_s / MACH)

    #dens, velx, vely, velz, pres = box.get_prims()
    dens, velx, vely, velz, pres = box.get_cons()

    fv = box.domainsize / np.array(dens.shape) # finite volume dimensions

    mach = np.sqrt(velx**2+vely**2+velz**2)/np.sqrt(pres/dens)
    vort = np.mean(fv)**5/12.0 * dens * ulz.norm(*ulz.curl(velx,vely,velz,fv[0],fv[1],fv[2]))

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
        data['time'], taskID,

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

def mkplot(taskID, ntasks, srcfp, sinkfp, crange=None):

    data = calc_data(srcfp)

    time = data['time']
    turntime = data['turntime']
    subplt = [2,2,0]
    fig = plt.figure(figsize=(20, 18))

    st = plt.suptitle(
        "decayturb periodic box: DG | t_d = % 2.4f (frame: %03d/%03d)" % (turntime, taskID+1, ntasks),
        fontsize='x-large')
    st.set_y(1.01)

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

    outfile = sinkfp % taskID
    plt.savefig(outfile,bbox_inches='tight')
    plt.close()

    print(outfile, flush=True)

with dslopts.Manager(scope=globals(), appendix="flashfiles can be defined after '--' or passed via stdin.") as mgr:
    mgr.add('sinkfp',  'path to store: <dir>/%03d.png')

srcfiles = list(_ignored_)

def task_minmax(x):
    taskID, srcfp = x
    return min_max(taskID, srcfp)

tmp = np.array(mpr.Pool().map(task_minmax,enumerate(srcfiles)))

i = ulz.mkincr(start=2)

crange = {
    'cmach': (np.min(tmp[:,next(i)]), np.max(tmp[:,next(i)])),
    'cdens': (np.min(tmp[:,next(i)]), np.max(tmp[:,next(i)])),
    'cpres': (np.min(tmp[:,next(i)]), np.max(tmp[:,next(i)])),
    'cvort': (np.min(tmp[:,next(i)]), np.max(tmp[:,next(i)]))
}

def task_plot(x): 
    taskID, srcfp = x
    return mkplot(taskID, len(srcfiles), srcfp, sinkfp, crange)

mpr.Pool().map(task_plot,enumerate(srcfiles))
