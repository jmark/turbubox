#!/usr/bin/env python

import numpy as np
import ulz
import interpolate
import flexi
import hopr

# func3d  = lambda x,y,z: np.sin(2*np.pi*x)
# func3d  = lambda x,y,z: np.ones_like(x)
# #func3d  = lambda x,y,z: x
# 
# lin1d = np.linspace(0.1,0.9,3)
# 
# xs = lin1d
# ys = lin1d
# zs = lin1d
# fs = func3d(*np.meshgrid(xs,ys,zs))
# 
# lin1dv = np.linspace(0.1,0.9,10)
# Xs,Ys,Zs = np.meshgrid(lin1dv,lin1dv,lin1dv,indexing='ij')
# 
# Fs = interpolate.lagrange_interpolate_3d(xs,ys,zs,fs,Xs,Ys,Zs).reshape(10,10,10)
# 
# print(np.all(np.abs(Fs-1.0) < 1e-8))

# xs = np.linspace(0,1,4)
# ms = np.meshgrid(*[xs]*2,indexing='ij')
# 
# Xs = np.linspace(0,1,10)
# Ms = np.meshgrid(*[Xs]*2,indexing='ij')
# 
# fs = np.ones_like(ms[0])
# 
# Fs = interpolate.lagrange_interpolate_2d_RG(xs,Xs,fs)

# xs = np.linspace(0,1,4)
# ms = np.meshgrid(*[xs]*3,indexing='ij')
# 
# Xs = np.linspace(0,1,10)
# Ms = np.meshgrid(*[Xs]*3,indexing='ij')
# 
# fs = np.ones_like(ms[0])
# 
# Fs = interpolate.lagrange_interpolate_3d_RG(xs,Xs,fs)

# flshfilepath = "/home/jmark/projects/stirturb/flexi-sims/whirlpool/turb-unit-128/flash_hdf5_chk_0050.h5"
# flexfilepath = "/home/jmark/projects/stirturb/flexi-sims/whirlpool/turb-unit-128/sim_State_0000000.000000000.h5"
# hoprfilepath = "/home/jmark/projects/stirturb/flexi-sims/whirlpool/turb-unit-128/sim_mesh.h5"


# flexfilepath = "/home/jmark/projects/stirturb/flexi-sims/whirlpool/constant/sim_State_0000000.000000000.h5"
# hoprfilepath = "/home/jmark/projects/stirturb/flexi-sims/whirlpool/constant/sim_mesh.h5"

flexfilepath = "/home/jmark/projects/stirturb/flexi-sims/whirlpool/plane/sim_State_0000000.000000000.h5"
hoprfilepath = "/home/jmark/projects/stirturb/flexi-sims/whirlpool/plane/sim_mesh.h5"



flx = flexi.File(flexfilepath, hopr.CartesianMeshFile(hoprfilepath))

xs = np.linspace(0,1,4)
Xs = np.linspace(0,1,5)

data = flx.data[:,:,:,:,0]

Fs = interpolate.flexi_to_box(xs, Xs, data, flx)

print(Fs.shape)
