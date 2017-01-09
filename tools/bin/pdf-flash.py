#!/usr/bin/env python3

# stdlib
import os, sys, pickle
import numpy as np

# jmark
import flash, ulz, dslopts
from defer_signals import DeferSignals

with dslopts.Manager(scope=globals(),appendix="flashfiles are be defined after '--'.") as mgr:
    mgr.add('sinkfptmpl', 'path template to store the pickle files: <dir>/03d%.pickle', str, '')
    mgr.add('usemultiproc', 'enable multiprocessing', dslopts.bool, True)
    mgr.add('skipfiles', 'skip already existing files', dslopts.bool, False)

def log(msg):
    print(msg, file=sys.stderr)

def pdf(data, nbins=1000):
    return np.histogram(data, bins=nbins, density=True)

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
    fls = flash.File(srcfp, 'r')

    # ndarrays
    dens = fls.data('dens')
    pres = fls.data('pres')
    velx = fls.data('velx')
    vely = fls.data('vely')
    velz = fls.data('velz')
    vels = np.sqrt(velx**2+vely**2+velz**2)
    ekin = np.prod(fls.cellsize)/2 * dens * (velx**2+vely**2+velz**2)
    vort = ulz.curl(velx,vely,velz,*tuple(fls.cellsize))

    # scalars
    time = fls.realscalars['time']
    step = fls.integerscalars['nstep']
    # snds = fls.realruntime['c_ambient']
    # dns0 = fls.realruntime['rho_ambient']
    # vrms = np.sqrt(np.sum(dens * (velx**2+vely**2+velz**2))/np.sum(dens))
    # mach = vrms / snds
    # dyns = snds * mach # dynamic speed
    # dynt = fls.domainsize[0] / dyns # dynamic time scale
    # turn = time / dynt # turn time

    res = {}
    res['taskid'] = taskid
    res['time'] = time
    res['step'] = step
    res['dens'] = pdf(dens, nbins='auto') 
    res['velx'] = pdf(velx)
    res['vely'] = pdf(vely)
    res['velz'] = pdf(velz)
    res['vels'] = pdf(vels)
    res['pres'] = pdf(pres)
    res['ekin'] = pdf(ekin)
    res['vort'] = pdf(vort)
    
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
