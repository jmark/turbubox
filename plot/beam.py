#!/usr/bin/env python3

import sys
from matplotlib import pylab as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import imp
import gausslobatto
import flash, flexi, hopr
import flash_to_flexi as flfl
import scipy.misc
import ulz
import interpolate
import glob

#PROJPATH = "/home/jmark/projects/stirturb/flexi-sims/smashing-balls"
#flshfilepath = "%s/flash_hdf5_chk_0080" % PROJPATH
#flexfilepath = "%s/run/sim_ERROR_State_0000000.787411102.h5" % PROJPATH

sys.argv.reverse()
progname = sys.argv.pop()
hoprfilepath = sys.argv.pop()
flexfilepath = sys.argv.pop()
sinkpath = sys.argv.pop() 

flx = flexi.File(flexfilepath, hopr.CartesianMeshFile(hoprfilepath))
cons = [flx.flexi_to_box(i) for i in range(0,8)]
prims = ulz.mhd_conservative_to_primitive(cons)
flx.h5file.close()

fig = plt.figure(figsize=(25, 12))

subplt = [2,4,0]

subplt[2] += 1
crange = {'vmin': -1, 'vmax': 1}
ys = np.sum(prims[1],axis=2).T
ys /= len(ys)
ax = fig.add_subplot(*subplt)
ax.set_title('FLEXI: column velx: %d^3' % len(ys))
ax.set_xlabel('x')
ax.set_ylabel('y')
#img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'),**crange)
img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'))
plt.colorbar(img,fraction=0.046, pad=0.04, format='%1.2f')

subplt[2] += 1
crange = {'vmin': 0, 'vmax': 2}
ys = np.sum(cons[4],axis=2).T
ys /= len(ys)
ax = fig.add_subplot(*subplt)
ax.set_title('FLEXI: column energy: %d^3' % len(ys))
ax.set_xlabel('x')
ax.set_ylabel('y')
#img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'),**crange)
img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'))
plt.colorbar(img,fraction=0.046, pad=0.04, format='%1.2f')

subplt[2] += 1
crange = {'vmin': 0, 'vmax': 2}
ys = np.sum(prims[4],axis=2).T
ys /= len(ys)
ax = fig.add_subplot(*subplt)
ax.set_title('FLEXI: column pressure: %d^3' % len(ys))
ax.set_xlabel('x')
ax.set_ylabel('y')
#img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'),**crange)
img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'))
plt.colorbar(img,fraction=0.046, pad=0.04, format='%1.2f')

subplt[2] += 1
crange = {'vmin': 0, 'vmax': 2}
ys = np.sum(prims[0],axis=2).T
ys /= len(ys)
ax = fig.add_subplot(*subplt)
ax.set_title('FLEXI: column density: %d^3' % len(ys))
ax.set_xlabel('x')
ax.set_ylabel('y')
#img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'),**crange)
img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'))
plt.colorbar(img,fraction=0.046, pad=0.04, format='%1.2f')


subplt[2] += 1
crange = {'vmin': -1, 'vmax': 1}
ys = np.sum(prims[1],axis=0).T
ys /= len(ys)
ax = fig.add_subplot(*subplt)
ax.set_title('FLEXI: column velx: %d^3' % len(ys))
ax.set_xlabel('x')
ax.set_ylabel('y')
#img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'),**crange)
img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'))
plt.colorbar(img,fraction=0.046, pad=0.04, format='%1.2f')

subplt[2] += 1
crange = {'vmin': 0, 'vmax': 2}
ys = np.sum(cons[4],axis=0).T
ys /= len(ys)
ax = fig.add_subplot(*subplt)
ax.set_title('FLEXI: column energy: %d^3' % len(ys))
ax.set_xlabel('x')
ax.set_ylabel('y')
#img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'),**crange)
img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'))
plt.colorbar(img,fraction=0.046, pad=0.04, format='%1.2f')

subplt[2] += 1
crange = {'vmin': 0, 'vmax': 2}
ys = np.sum(prims[4],axis=0).T
ys /= len(ys)
ax = fig.add_subplot(*subplt)
ax.set_title('FLEXI: column pressure: %d^3' % len(ys))
ax.set_xlabel('x')
ax.set_ylabel('y')
#img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'),**crange)
img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'))
plt.colorbar(img,fraction=0.046, pad=0.04, format='%1.2f')

subplt[2] += 1
crange = {'vmin': 0, 'vmax': 2}
ys = np.sum(prims[0],axis=0).T
ys /= len(ys)
ax = fig.add_subplot(*subplt)
ax.set_title('FLEXI: column density: %d^3' % len(ys))
ax.set_xlabel('x')
ax.set_ylabel('y')
#img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'),**crange)
img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'))
plt.colorbar(img,fraction=0.046, pad=0.04, format='%1.2f')



plt.savefig(sinkpath,bbox_inches='tight')
