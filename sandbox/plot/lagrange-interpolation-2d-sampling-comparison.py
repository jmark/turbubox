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

rms = lambda Fs: np.sqrt(np.sum(Fs**2)/len(Fs.ravel()))
rmse = lambda Fs,fs: np.sqrt(np.sum((Fs-fs)**2)/len(Fs.ravel()))
estr = lambda rms,rmse: "rms = %.2f | rms error = %.2f | relative rms error = %.2f%%" % (rms,rmse,rmse/rms * 100)


sys.argv.reverse()
sys.argv.pop()
outpath = sys.argv.pop()

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
visualpoints = np.meshgrid(*[np.linspace(-1,1,nvisu)]*ndim)   

def interpolate(nodes):
    # reference space
    samplepoints = np.meshgrid(*[nodes]*ndim)

    # mk interpolator in reference space
    Xs  = np.transpose(visualpoints).reshape(-1,ndim)
    ipl = gausslobatto.mk_lagrange_interpolator_2d(nodes,nodes,Xs)

    # interpolate in physical space
    ms = dom[0] + (dom[1]-dom[0]) * (np.array(samplepoints)+1)/2
    Ms = dom[0] + (dom[1]-dom[0]) * (np.array(visualpoints)+1)/2
    return ipl(f(*ms)).reshape(Ms[0].shape)


Ms = dom[0] + (dom[1]-dom[0]) * (np.array(visualpoints)+1)/2
Xs = np.transpose(Ms).reshape(-1,ndim)

Fs = f(*Ms)

fig = plt.figure(figsize=(18, 10))
ax = fig.add_subplot(2,2,1, projection='3d')
ax.set_title("original\n%s" % t,y=1.08)
if zlim: ax.set_zlim(*zlim)
ax.plot_wireframe(*Ms, Fs)

nodes = ulz.mk_body_centered_linspace(-1,1,npoly+1)
fs = interpolate(nodes)
e = estr(rms(Fs),rmse(Fs,fs))
ax = fig.add_subplot(2,2,2, projection='3d')
ax.set_title("body-centered sampling\n%s" % e,y=1.08)
if zlim: ax.set_zlim(*zlim)
ax.plot_wireframe(*Ms, fs)

nodes,_ = gausslobatto.LegendreGaussNodesAndWeights(npoly)
fs = interpolate(nodes)
e = estr(rms(Fs),rmse(Fs,fs))
ax = fig.add_subplot(2,2,3, projection='3d')
ax.set_title("Gauss nodes sampling\n%s" % e,y=1.08)
if zlim: ax.set_zlim(*zlim)
ax.plot_wireframe(*Ms, fs)

nodes,_ = gausslobatto.LegendreGaussLobattoNodesAndWeights(npoly)
fs = interpolate(nodes)
e = estr(rms(Fs),rmse(Fs,fs))
ax = fig.add_subplot(2,2,4, projection='3d')
ax.set_title("Gauss-Lobatto sampling\n%s" % e,y=1.08)
if zlim: ax.set_zlim(*zlim)
ax.plot_wireframe(*Ms, fs)

plt.savefig(outpath,bbox_inches='tight')
