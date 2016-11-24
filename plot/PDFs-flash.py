#!/usr/bin/env python3

# stdlib
import sys
import numpy as np
from itertools import chain
import multiprocessing
import pickle

# jmark
import flash, ulz, dslopts

def pdf(data):
    return np.histogram(data, bins=100, density=True)

def mkPDF(flashfp, sinkfp, taskID, ntasks=None):
    fls = flash.File(flashfp, mode='r')

    time = fls.realscalars['time']
    step = fls.integerscalars['nstep']
    c_s  = fls.realruntime['c_ambient']
    rho0 = fls.realruntime['rho_ambient']
    dynt = time / (fls.domainsize[0] / c_s / 10) # dynamic time

    GS = fls.gridsize
    DS = fls.domainsize
    CS = fls.cellsize

    Vgrid   = np.prod(GS)
    Vcell   = np.prod(CS) 
    Vdomain = np.prod(DS) 

    dens = fls.data('dens')
    pres = fls.data('pres')
    # velx = fls.data('velx')
    # vely = fls.data('velx')
    # velz = fls.data('velx')

    #mach = np.sqrt(ulz.norm(*vels)/3) / c_s
    #vort = CS[0]**5/12.0 * dens * ulz.norm(*ulz.curl(vels[0],vels[1],vels[2],CS[0],CS[1],CS[2]))

    result = {
        'time': time,
        'step': step,
        'c_s': c_s,
        'dens0': rho0,
        'dynt': dynt,
        'dens': pdf(dens)
    }

    outfile = sinkfp % taskID
    with open(outfile, 'wb') as fd:
        pickle.dump(result, fd)
    print(outfile)

with dslopts.Manager(scope=globals(),appendix="flashfiles can be defined after '--' or passed via stdin.") as mgr:
    mgr.add('sinkfilepath',  'path to store: <dir>/%03d.pickle')
    mgr.add('readfromstdin', 'yes/no', default='no')

if readfromstdin == 'yes':
    srcfiles = list(chain(_ignored_, map(str.strip,sys.stdin)))
else:
    srcfiles = list(_ignored_)

def task(x):
    mkPDF(x[1], sinkfilepath, x[0], len(srcfiles))

multiprocessing.Pool().map(task,enumerate(srcfiles))
