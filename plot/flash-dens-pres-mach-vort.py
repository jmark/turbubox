#!/usr/bin/env python3

# stdlib
import os
import sys
import numpy as np
from itertools import chain
import multiprocessing
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

matplotlib.rcParams.update({'font.size': 20})

# jmark
import flash, ulz, dslopts

def min_max(flashfp, taskID):
    fls  = flash.File(flashfp, mode='r')

    time = fls.realscalars['time']
    step = fls.integerscalars['nstep']
    c_s  = fls.realruntime['sona0']
    #rho0 = fls.realruntime['rho_ambient']
    LEN  = fls.domainsize[0]

    turntime = time / (LEN / c_s / 10)

    dens, velx, vely, velz, pres = fls.get_prims()
    mach = np.sqrt(velx**2+vely**2+velz**2)/np.sqrt(pres/dens)
    #ekin = fls.cellvolume/2 * np.sum(dens * (velx**2+vely**2+velz**2)) 

    CS = fls.cellsize
    vort = CS[0]**5/12.0 * dens * ulz.norm(*ulz.curl(velx,vely,velz,CS[0],CS[1],CS[2]))

    ax = 2
    cmach = np.log10(np.mean(mach,axis=ax))
    cdens = np.log10(np.mean(dens,axis=ax))
    cpres = np.log10(np.mean(pres,axis=ax))
    cvort = np.log10(np.mean(vort,axis=ax))

    min = np.min
    max = np.max

    return [
        time, taskID,

        min(cmach),
        max(cmach),

        min(cvort),
        max(cvort),

        min(cdens),
        max(cdens),

        min(cpres),
        max(cpres)
    ]

SOLVER = 'B5'
MACH = 10

def mkplot(flashfp, sinkfp, taskID, ntasks):

    outfile = sinkfp % taskID
    fls  = flash.File(flashfp, mode='r')

    time = fls.realscalars['time']
    step = fls.integerscalars['nstep']
    c_s  = fls.realruntime['sona0']
    #rho0 = fls.realruntime['rho_ambient']
    LEN  = fls.domainsize[0]

    turntime = time / (LEN / c_s / 10)

    dens, velx, vely, velz, pres = fls.get_prims()
    mach = np.sqrt(velx**2+vely**2+velz**2)/np.sqrt(pres/dens)
    #ekin = fls.cellvolume/2 * np.sum(dens * (velx**2+vely**2+velz**2)) 

    CS = fls.cellsize
    vort = CS[0]**5/12.0 * dens * ulz.norm(*ulz.curl(velx,vely,velz,CS[0],CS[1],CS[2]))

    ax = 2
    cdens = np.log10(np.sum(dens,axis=ax))
    cpres = np.log10(np.sum(pres,axis=ax))
    cvort = np.log10(np.sum(vort,axis=ax))
    cmach = np.mean(mach,axis=ax)

    subplt = [2,2,0]
    fig = plt.figure(figsize=(20, 18))

    st = plt.suptitle(
        "decay turb in periodic box: mach %d | %s | t_d = % 2.4f (frame: %03d/%03d)" % \
            (MACH, SOLVER, turntime, taskID+1, ntasks),
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

with dslopts.Manager(scope=globals(),appendix="flashfiles can be defined after '--' or passed via stdin.") as mgr:
    mgr.add('sinkfilepath',  'path to store: <dir>/%03d.png')
    mgr.add('solvername',    'name of solver: B3,B5,ES,...')
    mgr.add('machnumber',    '--', int)
    mgr.add('readfromstdin', 'yes/no', default='no')

srcfiles = []

if readfromstdin == 'yes':
    srcfiles = list(chain(_ignored_, map(str.strip,sys.stdin)))
else:
    srcfiles = list(_ignored_)

#srcfiles = list(filter(lambda i: os.path.isfile(sinkfilepath % i), range(len(flsfps))))

SOLVER = solvername
MACH   = machnumber

def task(x):
    return min_max(x[1], x[0])
    #return mkplot(x[1], sinkfilepath, x[0], len(srcfiles))

pool   = multiprocessing.Pool()
result = pool.map(task,enumerate(srcfiles))

for res in result:
    print(*res)
