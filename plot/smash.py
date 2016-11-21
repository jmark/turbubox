#!/usr/bin/env python3

# stdlib
import sys
import numpy as np
from itertools import chain
import multiprocessing
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 20})
from matplotlib import pyplot as plt

# jmark
import flexi, hopr, ulz, dslopts

title = \
'''Kelvin-Helmholtz Instability

  - FLEXI: 3rd order -> 4 nodes (GAUSS-LOBATTO) per element
  - periodic box: 64 x 64 x 1 elements
  - Euler equation, ideal monoatomic gas, isothermal EOS
  - subsonic: (relative) mach = 0.4 initial condition

frame: %03d/%03d'''

def trafo(data):
    axis = 2
    return np.sum(data,axis=axis).T/data.shape[axis]

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

def mkplot(meshfilepath, flexfilepath, sinkpath, taskID, ntasks):
    with flexi.File(flexfilepath, hopr.CartesianMeshFile(meshfilepath), mode='r') as flx:
        #TODO: distinguish navier and mhd case
        #cons = [flx.flexi_to_box(i) for i in range(0,8)]
        #prim = ulz.mhd_conservative_to_primitive(cons)
        cons = [flx.flexi_to_box(i) for i in range(0,5)]
        prim = ulz.navier_conservative_to_primitive(cons)

    subplt = [2,3,0]
    fig = plt.figure(figsize=(40, 22))

    def plot(data, title, crange=None):
        crange = None
        subplt[2] += 1
        ax = fig.add_subplot(*subplt)
        ax.set_title(title); ax.set_xlabel('x'); ax.set_ylabel('y')
        xs = np.arange(0,len(data)+4, 4)
        ax.set_xticks(xs); ax.set_yticks(xs); ax.grid()

        x0 = y0 = 0 - 1/2/len(data); x1 = y1 = len(data) + 1/2/len(data)
        args = {'X': data, 'cmap': plt.get_cmap('cubehelix'), 'extent': (x0,x1,y0,y1)}
        if crange is not None:
            args.update({'vmin': crange[0], 'vmax': crange[1]})
        img = ax.imshow(**args)

        plt.colorbar(img,fraction=0.045, pad=0.04, format='%1.2f')

    plot(trafo(prim[0]) ,'column density'  ,( 50,90))
    plot(trafo(prim[4]) ,'column pressure' ,( 45,100))
    plot(trafo(cons[4]) ,'column energy'   ,( 70,170))
    plot(trafo(prim[1]) ,'column velx'     ,(-15,32))
    plot(trafo(prim[2]) ,'column vely'     ,( 45,100))

    subplt[2] += 1
    ax = fig.add_subplot(*subplt)
    ax.axis('off')

    ax.text(0.0, 1.0,title % (taskID+1, ntasks),
        #fontsize='large',
        horizontalalignment='left',
        verticalalignment='top',
        transform = ax.transAxes)

    outfile = sinkpath % taskID
    plt.savefig(outfile,bbox_inches='tight')
    plt.close()
    print(outfile)

with dslopts.Manager(scope=globals(),appendix="flexifiles can be defined after '--' or passed via stdin.") as mgr:
    mgr.add('meshfilepath')
    mgr.add('sinkfilepath',  'path to store: <dir>/%03d.png')
    mgr.add('readfromstdin', 'yes/no', default='no')

if readfromstdin == 'yes':
    flsfps = chain(_ignored_, map(str.strip,sys.stdin))
else:
    flsfps = _ignored_

fps  = list(flsfps)
nfps = len(flsfps)

def task(x):
    #property(x[1], x[0])
    mkplot(meshfilepath, x[1], sinkfilepath, x[0], nfps)

multiprocessing.Pool().map(task,enumerate(fps))
