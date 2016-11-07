#!/usr/bin/env python

import numpy as np
import ulz
import interpolate
import flexi, hopr

flshfilepath = "/home/jmark/projects/stirturb/flexi-sims/whirlpool/turb-unit-128/flash_hdf5_chk_0050.h5"
flexfilepath = "/home/jmark/projects/stirturb/flexi-sims/whirlpool/turb-unit-128/sim_State_0000000.000000000.h5"
hoprfilepath = "/home/jmark/projects/stirturb/flexi-sims/whirlpool/turb-unit-128/sim_mesh.h5"

flx = flexi.File(flexfilepath, hopr.CartesianMeshFile(hoprfilepath))
fls = flash.File(flshfilepath)

box = np.zero(flx.mesh.gridsize)

npoly = flx.npoly
ntype = flx.nodetype
nvisu = npoly + 1

xs = np.linspace(0,1,nvisu)

Fs = interpolate.lagrange_interpolate_3d_RG(xs,Xs,fs)
