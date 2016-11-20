#!/usr/bin/env python3

# stdlib
import sys
import numpy as np
from itertools import chain
import multiprocessing
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

# jmark
import flexi, hopr, ulz, dslopts

def trafo(data):
    return np.sum(data,axis=2).T

def property(meshfile, flexfile, taskID):
    flx = flexi.File(flexfile, hopr.CartesianMeshFile(meshfile))
    cons = [flx.flexi_to_box(i) for i in range(0,8)]
    prim = ulz.mhd_conservative_to_primitive(cons)

    #func = lambda x: np.max(trafo(x))
    func = lambda x: np.min(trafo(x))
    print(
        func(prim[0]),
        func(prim[1]),
        func(prim[4]),
        func(cons[4]), flush=True)

def mkplot(meshfilepath, flexfilepath, sinkpath, taskID):
    with flexi.File(flexfilepath, hopr.CartesianMeshFile(meshfilepath), mode='r') as flx:
        cons = [flx.flexi_to_box(i) for i in range(0,8)]
        prim = ulz.mhd_conservative_to_primitive(cons)

    fig = plt.figure(figsize=(15, 12))
    subplt = [2,2,0]

    def plot(data, title, crange=None):
        subplt[2] += 1
        ax = fig.add_subplot(*subplt)
        ax.set_title('column %s: %d^3' % (title,len(data)))
        ax.set_xlabel('x'); ax.set_ylabel('y')
        if crange is not None:
            img = ax.imshow(data, cmap=plt.get_cmap('cubehelix'), vmin=crange[0], vmax=crange[1])
        else:
            img = ax.imshow(data, cmap=plt.get_cmap('cubehelix'))
        plt.colorbar(img,fraction=0.046, pad=0.04, format='%1.2f')

    plot(trafo(prim[0]), 'density', ( 50,90))
    plot(trafo(prim[1]), 'velx',    (-15,32))
    plot(trafo(prim[4]), 'pressure',( 45,100))
    plot(trafo(cons[4]), 'energy',  ( 70,170))

    outfile = sinkpath % taskID
    print(outfile)
    plt.savefig(outfile,bbox_inches='tight')
    plt.close()

with dslopts.Manager(scope=globals(),appendix="flexifiles can be defined after '--' or passed via stdin.") as mgr:
    mgr.add('meshfilepath')
    mgr.add('sinkfilepath',  'path to store: <dir>/%03d.png')
    mgr.add('readfromstdin', 'yes/no', default='no')

def task(x):
    #property(meshfilepath, x[1], x[0])
    mkplot(meshfilepath, x[1], sinkfilepath, x[0])

if readfromstdin == 'yes':
    flxfps = chain(_ignored_, map(str.strip,sys.stdin))
else:
    flxfps = _ignored_

multiprocessing.Pool().map(task,enumerate(flxfps))
