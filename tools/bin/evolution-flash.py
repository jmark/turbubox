#!/usr/bin/env python3

# stdlib
import os
import sys
import pickle
import numpy as np

# jmark
import flash, ulz, dslopts
from defer_signals import DeferSignals

with dslopts.Manager(scope=globals(),appendix="flashfiles are be defined after '--'.") as mgr:
    mgr.add('sinkfptmpl', 'path template to store the pickle files: <dir>/03d%.pickle', str, '')
    mgr.add('usemultiproc', 'enable multiprocessing', dslopts.bool, True)
    mgr.add('skipfiles', 'skip already existing files', dslopts.bool, False)

def evol(taskid, srcfp):
    # prepare sink file path
    try:
        snkfp = sinkfptmpl % taskid
    except TypeError:
        snkfp = sinkfptmpl

    # skip already done files
    if skipfiles and os.path.isfile(snkfp):
        return ['%s skipped.' % snkfp]

    # open flash file
    fls = flash.File(srcfp, 'r')

    # ndarrays
    dens = fls.data('dens')
    pres = fls.data('pres')
    velx = fls.data('velx')
    vely = fls.data('vely')
    velz = fls.data('velz')
    ekin = np.prod(fls.cellsize)/2 * dens * (velx**2+vely**2+velz**2)
    vort = ulz.curl(velx,vely,velz,*tuple(fls.cellsize))

    # scalars
    time = fls.realscalars['time']
    step = fls.integerscalars['nstep']
    snds = fls.realruntime['c_ambient']
    dns0 = fls.realruntime['rho_ambient']
    vrms = np.sqrt(np.sum(dens * (velx**2+vely**2+velz**2))/np.sum(dens))
    mach = vrms / snds
    dyns = snds * mach # dynamic speed
    dynt = fls.domainsize[0] / dyns # dynamic time scale
    turn = time / dynt # turn time

    ekintot = np.sum(ekin)
    vorttot = np.mean(fls.cellsize)**5/12.0 * np.sum(dens * (vort[0]**2 + vort[1]**2 + vort[2]**2))

    result = [taskid, step, time, turn, mach, ekintot, vorttot]
    
    if snkfp:
        with DeferSignals(): # make sure write is atomic
            with open(snkfp, 'wb') as fd:
                pickle.dump(result, fd)

    return result

srcfiles = map(str.rstrip, ARGV_TAIL)

if usemultiproc:
    from multiprocessing import Pool
    def task(x):
        return evol(x[0],x[1])
    Pool().map(task,enumerate(srcfiles))
else:
    for taskid, fp in enumerate(srcfiles):
        result = evol(taskid, fp)
        print(' '.join(map(str,result)), file=sys.stderr)
