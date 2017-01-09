#!/usr/bin/env python3

# stdlib
import os, sys, pickle
import numpy as np

# jmark
import flash, ulz, dslopts
from defer_signals import DeferSignals

from matplotlib import pylab as plt

with dslopts.Manager(scope=globals(),appendix="flashfiles are be defined after '--'.") as mgr:
    mgr.add('sinkfptmpl', 'path template to store the pickle files: <dir>/03d%.pickle', str, '')
    mgr.add('usemultiproc', 'enable multiprocessing', dslopts.bool, True)
    mgr.add('skipfiles', 'skip already existing files', dslopts.bool, False)

def log(msg):
    print(msg, file=sys.stderr)

def pdf(data, nbins=10000):
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
    velx = fls.data('velx')
    vely = fls.data('vely')
    velz = fls.data('velz')

    # scalars
    time = fls.realscalars['time']
    step = fls.integerscalars['nstep']

    Q = ulz.Q(velx,vely,velz).ravel()
    R = ulz.R(velx,vely,velz).ravel()

    nbins = 200  # number of bins in each dimension
    binrange = [[-0.1, 0.1], [-0.2, 0.1]]
    pdf, xbins, ybins = np.histogram2d(R,Q, bins=nbins, range=binrange)
    pdf /= pdf.sum()
    #pdf = np.log10(pdf)

    # i0 = (np.abs(xbins + 1)).argmin()
    # i1 = (np.abs(xbins - 1)).argmin()

    # j0 = (np.abs(ybins + 1)).argmin()
    # j1 = (np.abs(ybins - 1)).argmin()

    #print(i0,i1,j0,j1)

    #print(xbins.min(), xbins.max())
    #print(ybins.min(), ybins.max())

    #pdf = pdf[i0:i1,j0:j1]

    plt.figure(figsize=(10,10))
    #plt.imshow(pdf, origin='lower right')
    plt.imshow(pdf)
    plt.show()

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
