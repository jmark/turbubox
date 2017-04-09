#!/usr/bin/env pyturbubox

# stdlib
import os, sys, pickle
import numpy as np
from numpy.fft import rfftn, fftshift
import pathlib as pl

# jmark
import periodicbox, ulz
from shellavg import shell_avg_3d
from defer_signals import DeferSignals

def PositiveInt(arg):
    x = int(arg)
    if x >= 0:
        return x
    else:
        raise ValueError("'%d' must be positive!" % x)

def log(msg):
    print(msg, file=sys.stderr)

## ========================================================================= ##
## process commandline arguments

import argparse

pp = argparse.ArgumentParser(description = 'Batch Produce Powerspectra')

pp.add_argument(
    '--destdir',
    help='path to store: <dir>/%%03d.pickle',
    type=pl.Path, required=True,
)

pp.add_argument(
    '--parallel',
    help='enable parallel processes: -1 --> serial, 0 --> max. procs, N > 0 --> set n procs',
    type=int,
    default=-1,
)

# pp.add_argument(
#     '--normalize',
#     help='normalize',
#     action='store_true',
# )

pp.add_argument(
    '--skip',
    help='skip already produced files',
    action='store_true',
)

pp.add_argument(
    '--nsamples',
    help='number samples taken by shellavg3d',
    type=PositiveInt,
)

pp.add_argument(
    'snapshots',
    help='list of snapshot files',
    type=pl.Path,nargs='*',
)

ARGV = pp.parse_args()

## ========================================================================= ##
## Setup environment

if not ARGV.destdir.exists():
    ARGV.destdir.mkdir()

## ========================================================================= ##
## Routines

def mkpowerspectrum(taskid, srcfp):
    # open flash file
    fls = periodicbox.File(srcfp, mode='r')
    dens, velx, vely, velz, pres = fls.get_prims()

    vels = np.sqrt(velx**2+vely**2+velz**2)
    rhovels = dens**(1/3) * vels
    ekin = dens/2 * (velx**2+vely**2+velz**2)
    ekin = ekin / np.mean(ekin)
    #vort = ulz.curl(velx,vely,velz,*tuple(fls.cellsize))

    # scalars
    time = fls.time
    #step = fls.integerscalars['nstep']
    #snds = fls.realruntime['c_ambient']
    #dns0 = fls.realruntime['rho_ambient']
    #vrms = np.sqrt(np.sum(dens * (velx**2+vely**2+velz**2))/np.sum(dens))
    #mach = vrms / snds
    #dyns = snds * mach # dynamic speed
    #dynt = fls.domainsize[0] / dyns # dynamic time scale
    #turn = time / dynt # turn time

    fdens = fftshift(np.abs(rfftn(dens)))
    fvels = fftshift(np.abs(rfftn(vels)))
    fpres = fftshift(np.abs(rfftn(pres)))
    fekin = fftshift(np.abs(rfftn(ekin)))
    frhovels = fftshift(np.abs(rfftn(rhovels)))

    return dict(
        taskid = taskid,
        time = time,
        dens = shell_avg_3d(fdens**2, ARGV.nsamples),
        vels = shell_avg_3d(fvels**2, ARGV.nsamples),
        pres = shell_avg_3d(fpres**2, ARGV.nsamples),
        ekin = shell_avg_3d(fekin**2, ARGV.nsamples),
        rhovels = shell_avg_3d(frhovels**2, ARGV.nsamples),
    )

## ========================================================================= ##
## prepare task

def task(taskid, srcfp):
    snkfp = ARGV.destdir / srcfp.with_suffix('.pickle').name

    if ARGV.skip and snkfp.exists() and snkfp.stat().st_mtime > srcfp.stat().st_mtime:
        return
    else:
        log('Processing: ' + str(snkfp))
        retval = mkpowerspectrum(taskid, srcfp)
        with DeferSignals(): # make sure write is atomic
            with open(str(snkfp), 'wb') as fd:
                pickle.dump(retval, fd)

## ========================================================================= ##
## do work

if ARGV.parallel >= 0:
    import multiprocessing as mpr
    def _task(args):
        return task(*args)
    nprocs = None if ARGV.parallel == 0 else ARGV.parallel
    mpr.Pool(nprocs).map(_task,enumerate(ARGV.snapshots))
else:
    for i,x in enumerate(ARGV.snapshots):
        task(i,x)
