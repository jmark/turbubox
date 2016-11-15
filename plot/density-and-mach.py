#!/usr/bin/env python3

import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import ulz
import flash

matplotlib.rcParams.update({'font.size': 10})

sys.argv.reverse()
progname = sys.argv.pop()
flashfpB3 = sys.argv.pop()
flashfpB5 = sys.argv.pop()
ifile  = sys.argv.pop()
nfile  = sys.argv.pop()
sinkfp = sys.argv.pop()

fls  = flash.File(flashfpB3)
dens = fls.data('dens')

time = fls.realscalars['time']
step = fls.integerscalars['nstep']

c_s  = fls.realruntime['c_ambient']
rho0 = fls.realruntime['rho_ambient']

LEN  = fls.domainsize[0]

turntime = time / (LEN / c_s / 10)
fig = plt.figure(figsize=(15, 12))

subplt = [1,2,0]
#crange = {'vmin': 0, 'vmax': 4}
crange = {'vmin': 0, 'vmax': 4000}

subplt[2] += 1
#ys = np.log10(np.sum(dens,axis=2)).T
ys = np.sum(dens,axis=2).T
ax = fig.add_subplot(*subplt)
ax.set_title('column density: %s / t_c = %1.2f (frame: %s/%s)' % ('B3', turntime, ifile, nfile))
ax.set_xlabel('x index')
ax.set_ylabel('y index')
img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'),**crange)
#img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'))
plt.colorbar(img,fraction=0.046, pad=0.04, format='%1.2f')

fls  = flash.File(flashfpB5)
dens = fls.data('dens')

time = fls.realscalars['time']
step = fls.integerscalars['nstep']

c_s  = fls.realruntime['c_ambient']
rho0 = fls.realruntime['rho_ambient']

LEN  = fls.domainsize[0]

turntime = time / (LEN / c_s / 10)

subplt[2] += 1
#ys = np.log10(np.sum(dens,axis=2)).T
ys = np.sum(dens,axis=2).T
ax = fig.add_subplot(*subplt)
ax.set_title('column density: %s / t_c = %1.2f (frame: %s/%s)' % ('B5', turntime, ifile, nfile))
ax.set_xlabel('x index')
ax.set_ylabel('y index')
img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'),**crange)
#img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'))
plt.colorbar(img,fraction=0.046, pad=0.04, format='%1.2f')

plt.savefig(sinkfp,bbox_inches='tight')