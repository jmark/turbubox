#!/usr/bin/env python3

import numpy as np
from numpy.fft import rfftn, fftshift
import sys
import pathlib

import flash
from shellavg import shell_avg_3d
import ulz
import dslopts

def ExistingPath(arg):
    pth = pathlib.Path(arg)
    if not pth.exists():
        raise OSError("'%s' does not exists!" % pth)
    return pth

def PositiveInt(arg):
    x = int(arg)
    if x >= 0:
        return x
    raise ValueError("'%d' must be positive!" % x)

def evol(taskid, srcfp, snkfptmpl=''):

    snkfp = snkfptmpl % taskid
    fls   = flash.File(srcfp)

    # ndarrays
    dens = fls.data('dens')
    pres = fls.data('pres')
    velx = fls.data('velx')
    vely = fls.data('vely')
    velz = fls.data('velz')
    vels = np.sqrt(velx**2+vely**2+velz**2)
    rhovels = dens**(1/3) * vels
    ekin = np.prod(fls.cellsize)/2 * dens * (velx**2+vely**2+velz**2)
    #vort = ulz.curl(velx,vely,velz,*tuple(fls.cellsize))

    # scalars
    time = fls.realscalars['time']
    step = fls.integerscalars['nstep']
    snds = fls.realruntime['c_ambient']
    dns0 = fls.realruntime['rho_ambient']
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
    res['time'] = time
    res['step'] = step
    res['dens'] = shell_avg_3d(fdens**2, nsamples)
    res['vels'] = shell_avg_3d(fvels**2, nsamples)
    res['pres'] = shell_avg_3d(fpres**2, nsamples)
    res['ekin'] = shell_avg_3d(fekin**2, nsamples)
    res['rhovels'] = shell_avg_3d(frhovels**2, nsamples)
    
    if snkfptmpl:
        with open(snkfp, 'wb') as fd:
            pickle.dump(res, fd)

with dslopts.Manager(scope=globals(),appendix="flashfiles are be defined after '--'.") as mgr:
    mgr.add('sinkfp', 'path to store the pickle files: <dir>/03d%.pickle', str, '')
    mgr.add('nsamples' ,'no. of samples'  ,PositiveInt, 0)
    mgr.add('usemultiproc', 'enable multiprocessing', dslopts.bool, True)
    mgr.add('skipfiles', 'skip already existing files', dslopts.bool, True)

srcfiles = ARGVTAIL

if usemultiproc:
    from multiprocessing import Pool
    def task(x):
        return evol(x[0],x[1], sinkfp)
    Pool().map(task,enumerate(srcfiles))
else:
    for taskid, fp in enumerate(srcfiles):
        result = evol(taskid, fp, sinkfp)
        print(' '.join(map(str,result)), file=sys.stderr)
