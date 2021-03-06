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

def property(flashfp, taskID):
    fls = flash.File(flashfp, mode='r')

    time = fls.realscalars['time']
    step = fls.integerscalars['nstep']
    c_s  = fls.realruntime['c_ambient']
    rho0 = fls.realruntime['rho_ambient']
    LEN  = fls.domainsize[0]

    turntime = time / (LEN / c_s / 10)

    GS = fls.gridsize
    DS = fls.domainsize
    CS = fls.cellsize

    Vgrid   = np.prod(GS)
    Vcell   = np.prod(CS) 
    Vdomain = np.prod(DS) 

    dens = fls.data('dens')
    pres = fls.data('pres')
    vels = [fls.data('vel'+dim) for dim in 'x y z'.split()]
    mach = np.sqrt(ulz.norm(*vels)/3) / c_s
    vort = CS[0]**5/12.0 * dens * ulz.norm(*ulz.curl(vels[0],vels[1],vels[2],CS[0],CS[1],CS[2]))

    ax = 2
    cdens = np.log10(np.sum(dens,axis=ax))
    cpres = np.log10(np.sum(pres,axis=ax))
    cmach = np.sum(mach,axis=ax)/mach.shape[ax]
    cvort = np.log10(np.sum(vort,axis=ax))

    #func = np.min
    func = np.max
    print(
        func(cdens),
        func(cpres),
        func(cmach),
        func(cvort), flush=True)

SOLVER = 'B5'
MACH = 10

def mkplot(flashfp, sinkfp, taskID, ntasks):

    fls = flash.File(flashfp, mode='r')

    time = fls.realscalars['time']
    step = fls.integerscalars['nstep']
    c_s  = fls.realruntime['c_ambient']
    rho0 = fls.realruntime['rho_ambient']
    LEN  = fls.domainsize[0]

    turntime = time / (LEN / c_s / 10)

    GS = fls.gridsize
    DS = fls.domainsize
    CS = fls.cellsize

    Vgrid   = np.prod(GS)
    Vcell   = np.prod(CS) 
    Vdomain = np.prod(DS) 

    dens = fls.data('dens')
    pres = fls.data('pres')
    velx = fls.data('velx')
    vely = fls.data('vely')
    velz = fls.data('velz')
    mach = np.sqrt((velx**2+vely**2+velz**2)/3)

    ekin = Vcell/2 * np.sum(dens * (velx**2+vely**2+velz**2)) 

    vort = CS[0]**5/12.0 * dens * ulz.norm(*ulz.curl(velx,vely,velz,CS[0],CS[1],CS[2]))

    ax = 2
    cdens = np.log10(np.sum(dens,axis=ax))
    cpres = np.log10(np.sum(pres,axis=ax))
    cvelx = np.sum(velx,axis=ax)
    cvely = np.sum(vely,axis=ax)
    cmach = np.sum(mach,axis=ax)/mach.shape[ax]
    cvort = np.log10(np.sum(vort,axis=ax))

    subplt = [2,2,0]
    fig = plt.figure(figsize=(20, 18))

    st = plt.suptitle(
        "Stirred turbulence in periodic box: mach %d | %s | t_d = % 2.4f (frame: %03d/%03d)" % \
            (MACH, SOLVER, turntime, taskID+1, ntasks),
        fontsize='x-large')
    st.set_y(1.01)

    def plot(data, title, crange=None):
        subplt[2] += 1
        ax = fig.add_subplot(*subplt)
        ax.set_title(title)
        ax.set_xlabel('x index'); ax.set_ylabel('y index')
        if crange is not None:
            img = ax.imshow(data, cmap=plt.get_cmap('cubehelix'), vmin=crange[0], vmax=crange[1])
        else:
            img = ax.imshow(data, cmap=plt.get_cmap('cubehelix'))
        plt.colorbar(img,fraction=0.0456, pad=0.04, format='%1.2f')

    plot(cdens, 'column density (log10)', (0,4))
    plot(cpres, 'column pressure (log10)', (0,4))
    plot(cvelx, 'column velocity x', (-2,2))
    plot(cvely, 'column velocity y', (-2,2))
    #plot(cmach, 'column sonic mach number (grid normalized)', (0,10))
    #plot(cvort, 'column vorticity (log10)', (-10,-4))

    outfile = sinkfp % taskID
    fig.tight_layout()
    plt.savefig(outfile,bbox_inches='tight')
    plt.close()
    print(outfile)

with dslopts.Manager(scope=globals(),appendix="flashfiles can be defined after '--' or passed via stdin.") as mgr:
    mgr.add('sinkfilepath',  'path to store: <dir>/%03d.png')
    mgr.add('solvername',    'name of solver: B3,B5,ES,...')
    mgr.add('machnumber',    '--', int)
    mgr.add('readfromstdin', 'yes/no', default='no')

if readfromstdin == 'yes':
    flsfps = list(chain(_ignored_, map(str.strip,sys.stdin)))
else:
    flsfps = list(_ignored_)

srcfiles = list(filter(lambda i: os.path.isfile(sinkfilepath % i), range(len(flsfps))))

SOLVER = solvername
MACH   = machnumber

def task(x):
    #property(x[1], x[0])
    mkplot(x[1], sinkfilepath, x[0], len(srcfiles))

multiprocessing.Pool().map(task,enumerate(flsfps))
