#!/usr/bin/env python3

# stdlib
import os
import sys
import numpy as np
#import multiprocessing
#from itertools import chain

# import matplotlib
# matplotlib.use('Agg')
# from matplotlib import pyplot as plt
# 
# matplotlib.rcParams.update({'font.size': 20})

# jmark
import flexi, ulz, dslopts

SOLVER = 'DG'
MACH = 2

def min_max(srcfp, taskID):
    box = flexi.PeriodicBox(srcfp)

    time = box.time
    c_s  = 1
    turntime = time / (np.mean(box.domainsize) / c_s / MACH)

    dens, velx, vely, velz, pres = box.get_prims()

    fv = box.domainsize / np.array(dens.shape) # finite volume dimensions

    mach = np.sqrt(velx**2+vely**2+velz**2)/np.sqrt(pres/dens)
    vort = np.mean(fv)**5/12.0 * dens * ulz.norm(*ulz.curl(velx,vely,velz,fv[0],fv[1],fv[2]))
    #ekin = fv**3/2 * np.sum(dens * (velx**2+vely**2+velz**2)) 

    ax = 2
    cdens = np.log10(np.sum(dens,axis=ax))
    cpres = np.log10(np.sum(pres,axis=ax))
    cvort = np.log10(np.sum(vort,axis=ax))
    cmach = np.mean(mach,axis=ax)

    min = np.min
    max = np.max

    return [
        time, taskID,

        min(cdens),
        max(cdens),

        min(cpres),
        max(cpres),

        min(cvort),
        max(cvort),

        min(cmach),
        max(cmach)
    ]

def mkplot(srcfp, meshfp, sinkfp, taskID, ntasks):

    outfile = sinkfp % taskID
    box = flexi.PeriodicBox(srcfp, meshfp)

    time = box.time
    LEN  = box.hopr.domainsize[0]
    c_s  = 1
    turntime = time / (LEN / c_s / MACH)

    dens, velx, vely, velz, pres = box.get_prims()
    mach = np.sqrt(velx**2+vely**2+velz**2)/np.sqrt(pres/dens)
    #ekin = box.cellvolume/2 * np.sum(dens * (velx**2+vely**2+velz**2)) 

    CS = box.hopr.domainsize / np.array(dens.shape)
    vort = CS[0]**5/12.0 * dens * ulz.norm(*ulz.curl(velx,vely,velz,CS[0],CS[1],CS[2]))

    ax = 2
    cdens = np.log10(np.sum(dens,axis=ax))
    cpres = np.log10(np.sum(pres,axis=ax))
    cvort = np.log10(np.sum(vort,axis=ax))
    cmach = np.mean(mach,axis=ax)

    subplt = [2,2,0]
    fig = plt.figure(figsize=(20, 18))

    st = plt.suptitle(
        "decayturb in periodic box: mach %d | %s | t_d = % 2.4f (frame: %03d/%03d)" % \
            (MACH, SOLVER, 1/MACH, taskID+1, ntasks),
        fontsize='x-large')
    st.set_y(1.01)

    def plot(data, title, crange=None):
        subplt[2] += 1
        ax = fig.add_subplot(*subplt)
        ax.set_title(title)
        ax.set_xlabel('x index'); ax.set_ylabel('y index')
        if crange is not None:
            img = ax.imshow(data, cmap=plt.get_cmap('cubehelix'), interpolation='none', vmin=crange[0], vmax=crange[1])
        else:
            img = ax.imshow(data, cmap=plt.get_cmap('cubehelix'), interpolation='none')
        plt.colorbar(img,fraction=0.0456, pad=0.04, format='%1.2f')

    crange = None

    crange = (2.0,3.0)
    plot(cdens, 'column density (log10)', crange)

    crange = (2.0,3.0) 
    plot(cpres, 'column pressure (log10)', crange)

    crange = (-8.5,-6) 
    plot(cvort, 'column vorticity (log10)', crange)

    crange = (0.5,3) 
    plot(cmach, 'column sonic mach number (grid normalized)', crange)

    fig.tight_layout()
    plt.savefig(outfile,bbox_inches='tight')
    plt.close()
    print(outfile)

with dslopts.Manager(scope=globals(), appendix="flashfiles can be defined after '--' or passed via stdin.") as mgr:
    mgr.add('sinkfp',  'path to store: <dir>/%03d.png')
    #mgr.add('meshfp',  'path to mesh', str, '')

flsfps = list(_ignored_)
srcfiles = list(filter(lambda i: os.path.isfile(sinkfp % i), range(len(flsfps))))

def task(x):
    return min_max(x[1], x[0])
    #return mkplot(x[1], sinkfp, x[0], len(srcfiles))

#pool   = multiprocessing.Pool()
#result = pool.map(task,enumerate(srcfiles))

for fp in srcfiles:
    res = min_max(fp, 0) 
    print(*res)
