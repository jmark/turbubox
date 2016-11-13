#!/usr/bin/env python3

import sys
from matplotlib import pyplot as plt
#from mpl_toolkits.mplot3d import axes3d
import numpy as np
import ulz
import flash

sys.argv.reverse()
progname = sys.argv.pop()
flashfpB3 = sys.argv.pop()
flashfpB5 = sys.argv.pop()
sinkfp = sys.argv.pop()

flsB3 = flash.File(flashfpB3)
flsB5 = flash.File(flashfpB5)

densityB3   = flsB3.data('dens')
densityB5   = flsB5.data('dens')

fig = plt.figure(figsize=(15, 12))

subplt = [1,2,0]
crange = None

subplt[2] += 1
#crange = {'vmin': -1, 'vmax': 1}
ys = np.sum(densityB3,axis=2).T
ys /= len(ys)
ax = fig.add_subplot(*subplt)
ax.set_title('B3 - column density: %d^3' % len(ys))
ax.set_xlabel('x')
ax.set_ylabel('y')
#img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'),**crange)
img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'))
plt.colorbar(img,fraction=0.046, pad=0.04, format='%1.2f')

subplt[2] += 1
#crange = {'vmin': 0, 'vmax': 2}
ys = np.sum(densityB5,axis=2).T
ys /= len(ys)
ax = fig.add_subplot(*subplt)
ax.set_title('B5 - column density: %d^3' % len(ys))
ax.set_xlabel('x')
ax.set_ylabel('y')
#img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'),**crange)
img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'))
plt.colorbar(img,fraction=0.046, pad=0.04, format='%1.2f')

plt.savefig(sinkfp,bbox_inches='tight')
