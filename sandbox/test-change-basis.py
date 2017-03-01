#!/usr/bin/env python3

import sys
sys.path.extend(["/home/jmark/projects/stirturb/turbubox/tools/%s" % s for s in "bin lib".split()])

import numpy as np
import gausslobatto
import flash, flexi, hopr
import ulz
import interpolate

gamma = 5/3

srcfp = "/mnt/seagate/flexi/stirturb/run-shockcapturing/sim_State_0000000.350000000.h5"
#srcfp = "/home/jmark/dump/data/sim_State_0000002.450000000.h5"
#srcfp = "/mnt/seagate/flexi/stirturb/snapshots/sim_State_0000000.500000000.h5"

with flexi.PeriodicBox(srcfp) as ff:
    #dens,momx,momy,momz,ener = ff.get_cons_fv()
    dens,velx,vely,velz,pres = ff.get_prims_fv()
    npoly = ff.Npoly

xs = gausslobatto.mk_nodes_from_to(-1,1,npoly)
ys = xs

# modal Vandermonde Matrix
Vmodal = np.zeros((npoly+1,npoly+1))
for i in range(0,npoly+1):
    for j in range(0,npoly+1):
        Vmodal[i,j] = gausslobatto.LegendrePolynomialAndDerivative(j,xs[i])[0]
Wmodal = np.linalg.inv(Vmodal)

modal = interpolate.change_basis(Wmodal, pres)
