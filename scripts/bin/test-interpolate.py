#!/usr/bin/env python

import sys
sys.path.extend(["/home/jmark/data/projects/stirturb/turbubox/pyflash/%s" % s for s in "bin lib cython".split()])

import numpy as np
import ulz
#import gausslobatto
import interpolate

func3d  = lambda x,y,z: np.sin(2*np.pi*x)
func3d  = lambda x,y,z: np.ones_like(x)
#func3d  = lambda x,y,z: x

lin1d = np.linspace(0.1,0.9,3)

xs = lin1d
ys = lin1d
zs = lin1d
fs = func3d(*np.meshgrid(xs,ys,zs))

lin1dv = np.linspace(0.1,0.9,10)
Xs,Ys,Zs = np.meshgrid(lin1dv,lin1dv,lin1dv,indexing='ij')

Fs = interpolate.lagrange_interpolate_3d(xs,ys,zs,fs,Xs,Ys,Zs).reshape(10,10,10)

print(np.all(np.abs(Fs-1.0) < 1e-8))
