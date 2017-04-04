#!/usr/bin/env pyturbubox

# stdlib
import os, sys, pickle
import numpy as np
from numpy.fft import rfftn, fftshift

# jmark
import periodicbox, ulz, dslopts
from shellavg import shell_avg_3d
from defer_signals import DeferSignals

def PositiveInt(arg):
    x = int(arg)
    if x >= 0:
        return x
    raise ValueError("'%d' must be positive!" % x)

def log(msg):
    print(msg, file=sys.stderr)

with dslopts.Manager(scope=globals(),appendix="flashfiles are be defined after '--'.") as mgr:
    mgr.add('sinkfptmpl', 'path template to store the pickle files: <dir>/03d%.pickle', str, '')
    mgr.add('nsamples' ,'no. of samples'  ,PositiveInt, 0)
    mgr.add('usemultiproc', 'enable multiprocessing', int, 1)
    mgr.add('skipfiles', 'skip already existing files', int, 0)

def task(taskid, srcfp):
    # prepare sink file path
    try:
        snkfp = sinkfptmpl % taskid
    except TypeError:
        snkfp = sinkfptmpl

    # skip already done files
    if skipfiles and os.path.isfile(snkfp):
        log('%s skipped.' % snkfp)
        return

    # open flash file
    fls = periodicbox.File(srcfp, 'r')
    dens, velx, vely, velz, pres = fls.get_prims()

    vels = np.sqrt(velx**2+vely**2+velz**2)
    rhovels = dens**(1/3) * vels
    ekin = np.prod(fls.cellsize)/2 * dens * (velx**2+vely**2+velz**2)
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

    res = {}
    res['taskid'] = taskid
    res['time'] = time
    #res['step'] = step
    res['dens'] = shell_avg_3d(fdens**2, nsamples)
    res['vels'] = shell_avg_3d(fvels**2, nsamples)
    res['pres'] = shell_avg_3d(fpres**2, nsamples)
    res['ekin'] = shell_avg_3d(fekin**2, nsamples)
    res['rhovels'] = shell_avg_3d(frhovels**2, nsamples)
    
    if snkfp:
        with DeferSignals(): # make sure write is atomic
            with open(snkfp, 'wb') as fd:
                pickle.dump(res, fd)

    log(snkfp)

srcfiles = map(str.rstrip, ARGV_TAIL)

if usemultiproc:
    from multiprocessing import Pool
    def _task(x):
        return task(x[0],x[1])
    Pool().map(_task,enumerate(srcfiles))
else:
    for taskid, fp in enumerate(srcfiles):
        task(taskid, fp)
