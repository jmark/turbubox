#!/usr/bin/env python3

# stdlib
import sys
import pickle
import numpy as np
import multiprocessing

# jmark
import flash, ulz, dslopts

def task(taskid, fp):
    fls = flash.File(fp)

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

    result = [count, step, time, turn, mach, ekintot, vorttot]
    print("\t".join(map(str,result)), file=sys.stderr)
    return result

with dslopts.Manager(scope=globals(),appendix="flashfiles are be defined after '--'.") as mgr:
    mgr.add('sinkfilepath',  'path to store the pickle file')

srcfiles = _ignored_
result = multiprocessing.Pool().map(lambda x: task(*x),enumerate(srcfiles))

with open(sinkfilepath, 'wb') as fd:
    pickle.dump(result, fd)
