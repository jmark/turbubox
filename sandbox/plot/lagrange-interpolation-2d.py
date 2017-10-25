#!/usr/bin/env python3

import sys
from matplotlib import pylab as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import gausslobatto
import flash, flexi, hopr
import flash_to_flexi as flfl
import scipy.misc
import ulz

sys.argv.reverse()
sys.argv.pop()
outpath = sys.argv.pop()

f = lambda x,y: x+y
t = "f(x,y) = x + y"
zlim = None

f = lambda x,y: x**2 + x*y - y**2
t = "f(x,y) = x^2 - x.y - y^2"
zlim = None

f = lambda x,y: -np.sin(x/4)-np.cos(y/4)
t = "f(x,y) = -(sin(x/4)+cos(y/4))"
zlim = None

f = lambda x,y: -np.sin(x/2)-np.cos(y/4)
t = "f(x,y) = -(sin(x/2)+cos(y/4))"
zlim = (-1.5,1.5)

# f = lambda x,y: -(np.sin(x/2)+3*np.cos(y))
# t = "f(x,y) = -(100*sin(x/2)+cos(y/4))"
# zlim = (-4,4)

npoly = 3
ndim  = 2
nvisu = 20
dom   = (-2*np.pi,2*np.pi)

# reference space
nodes, weights = gausslobatto.LegendreGaussLobattoNodesAndWeights(npoly)
samplepoints = np.meshgrid(*(nodes,)*ndim)
visualpoints = np.meshgrid(*(np.linspace(-1,1,nvisu),)*ndim)   

# mk interpolator
Xs  = np.transpose(visualpoints).reshape(-1,ndim)
ipl = gausslobatto.mk_lagrange_interpolator_2d(nodes,nodes,Xs)

# transform from ref space to physical space
ms = dom[0] + (dom[1]-dom[0]) * (np.array(samplepoints)+1)/2
Ms = dom[0] + (dom[1]-dom[0]) * (np.array(visualpoints)+1)/2

# get original and interpolated values
fs = ipl(f(*ms)).reshape(Ms[0].shape)
Fs = f(*Ms)

# estimate error
rms = np.sqrt(np.sum(Fs**2)/len(Fs.ravel()))
rmse = np.sqrt(np.sum((Fs-fs)**2)/len(Fs.ravel()))
estring = "rms = %f | rms error = %f | relative rms error = %f%%" % (rms,rmse,rmse/rms * 100)
print(estring)

# plot
fig = plt.figure(figsize=(18, 5))

ax = fig.add_subplot(1,2,1, projection='3d')
ax.set_title("original: %s" % t,y=1.08)
if zlim: ax.set_zlim(*zlim)
ax.plot_wireframe(*Ms, Fs)

ax = fig.add_subplot(1,2,2, projection='3d')
ax.set_title("interpolated",y=1.08)
if zlim: ax.set_zlim(*zlim)
ax.plot_wireframe(*Ms, fs)

plt.savefig(outpath,bbox_inches='tight')
