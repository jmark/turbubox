#!/usr/bin/env python

import numpy as np
import ulz
import interpolate
import flexi
import hopr
import flash

if 0:
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

if 0:
    xs = np.linspace(0,1,4)
    ms = np.meshgrid(*[xs]*2,indexing='ij')

    Xs = np.linspace(0,1,10)
    Ms = np.meshgrid(*[Xs]*2,indexing='ij')

    fs = np.ones_like(ms[0])

    Fs = interpolate.lagrange_interpolate_2d_RG(xs,Xs,fs)

if 0:
    xs = np.linspace(0,1,4)
    ms = np.meshgrid(*[xs]*3,indexing='ij')

    Xs = np.linspace(0,1,10)
    Ms = np.meshgrid(*[Xs]*3,indexing='ij')

    fs = np.ones_like(ms[0])

    Fs = interpolate.lagrange_interpolate_3d_RG(xs,Xs,fs)

# flshfilepath = "/home/jmark/projects/stirturb/flexi-sims/whirlpool/turb-unit-128/flash_hdf5_chk_0050.h5"
# flexfilepath = "/home/jmark/projects/stirturb/flexi-sims/whirlpool/turb-unit-128/sim_State_0000000.000000000.h5"
# hoprfilepath = "/home/jmark/projects/stirturb/flexi-sims/whirlpool/turb-unit-128/sim_mesh.h5"

# flexfilepath = "/home/jmark/projects/stirturb/flexi-sims/whirlpool/constant/sim_State_0000000.000000000.h5"
# hoprfilepath = "/home/jmark/projects/stirturb/flexi-sims/whirlpool/constant/sim_mesh.h5"

if 0:
    flx = flexi.File(flexfilepath, hopr.CartesianMeshFile(hoprfilepath))

    xs = np.linspace(0,1,4)
    Xs = np.linspace(0,1,5)

    data = flx.data[:,:,:,:,0]

    Fs = interpolate.flexi_to_box(xs, Xs, data, flx)

    print(Fs.shape)

if 1:
    flexfilepath = "/home/jmark/projects/stirturb/flexi-sims/whirlpool/plane/sim_State_0000000.000000000.h5"
    hoprfilepath = "/home/jmark/projects/stirturb/flexi-sims/whirlpool/plane/sim_mesh.h5"

    flx = flexi.File(flexfilepath, hopr.CartesianMeshFile(hoprfilepath))
    fls = flash.FakeFile([[0,0,0],[1,1,1]],(flx.mesh.gridsize * (flx.npoly+1)), fillby='constant')

    xs = np.linspace(0,1,4)
    ms = np.meshgrid(*[xs]*3,indexing='ij')

    Xs = np.linspace(0,1,4)
    Ms = np.meshgrid(*[Xs]*3,indexing='ij')

    #box = ulz.wrap_in_guard_cells(fls.data('dens'))
    box = fls.data('dens')
    box = np.ones_like(box)
    #Fs = interpolate.box_to_flexi_with_averaged_boundaries(xs,Xs,box,flx)

    elems = interpolate.box_to_elements_avg_boundaries(box,flx)

    print(elems.shape)
    print(np.all(np.abs(box - 1.0) < 1e-8))
    print(np.all(np.abs(elems - 1.0) < 1e-8))
    #print(Fs[-1]) 

if 0:
    #arr = np.ones([3,3,3])
    arr = np.arange(0,4*5*6).reshape(4,5,6)
    interpolate.test(arr)
