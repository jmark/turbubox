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

case = 4

if case == 0:
    casename = 'plane'
    fillby   = 'planeX'
    stride   = 1
    zlim     = (0,1)

elif case == 1:
    casename = 'plane+wiggle'
    fillby   = 'plane+wiggle'
    stride   = 2
    zlim      = (0,3)

elif case == 2:
    casename = 'gaussian'
    fillby   = 'gaussianXYZ'
    stride   = 1
    zlim     = (0,1)

elif case == 3:
    casename = 'steps'
    fillby   = casename+'XYZ'
    stride   = 2
    zlim     = (0,5)

elif case == 4:
    casename = 'turb-unit-128'
    fillby   = 'file'
    stride   = 2
    zlim     = (0,5)



rms  = lambda Fs: np.sqrt(np.sum(Fs**2)/len(Fs.ravel()))
rmse = lambda Fs,fs: np.sqrt(np.sum((Fs-fs)**2)/len(Fs.ravel()))
estr = lambda rms,rmse: "rms = %f | rms error = %f | relative rms error = %f%%" % (rms,rmse,rmse/rms * 100)

sys.argv.reverse()
sys.argv.pop()
outpath = sys.argv.pop()

flexfilepath = "/home/jmark/data/projects/stirturb/flexi-sims/whirlpool/%s/sim_State_0000000.000000000.h5" % casename
hoprfilepath = "/home/jmark/data/projects/stirturb/flexi-sims/whirlpool/%s/sim_mesh.h5" % casename
flshfilepath = "/home/jmark/data/projects/stirturb/flexi-sims/whirlpool/%s/flash_hdf5_chk_0050.h5" % casename

flx = flexi.File(flexfilepath, hopr.CartesianMeshFile(hoprfilepath))

if 'file' in fillby:
    fls = flash.File(flshfilepath)
else:
    fls = flfl.FakeFlashFile(flx.mesh.domain,flx.mesh.gridsize.astype(int)*(flx.npoly+1),fillby=fillby)

flsdata = fls.data('dens')

## INTERPOLATION ##

npoly = flx.npoly
ntype = flx.nodetype
nvisu = 4

#RGcube = ulz.mk_cartesian_product_3d(*[ulz.mk_body_centered_linspace(-1,1,nvisu)]*3)
RGcube = ulz.mk_cartesian_product_3d(*[np.linspace(-1+0.0001,1-0.0001,nvisu)]*3)
RGgrid = np.array([ll + (tr-ll)*(RGcube+1)/2 for ll,tr in zip(*flx.mesh.get_cell_coords())]).reshape(-1, *[nvisu]*3, 3)
ipl    = gausslobatto.mk_lagrange_interpolator_3d(*[gausslobatto.mk_nodes(npoly, ntype)]*3, RGcube)

flxdata = np.array([ipl(f).reshape((nvisu,)*3) for f in flx.data[:,:,:,:,0].transpose(0,3,2,1)])
flxgrid, flxdata = ulz.sort_unstructured_grid(RGgrid, flxdata)

#print(estr(rms(flsdata),rmse(flsdata,flxdata)))

## PLOTTING ##

flslin1d = ulz.mk_body_centered_linspace(0,1, len(flsdata))
flxlin1d = flxgrid.T[0,0,0]  

plt.figure(figsize=(12,6))

plt.grid()
plt.xticks(np.linspace(0,1,flx.mesh.gridsize[0],endpoint=False))

for label in plt.gca().xaxis.get_ticklabels():
    label.set_visible(False)

for label in plt.gca().xaxis.get_ticklabels()[::3]:
    label.set_visible(True)

i = len(flsdata)//2
j = len(flsdata)//2 
k = len(flsdata)//2 

xs = flslin1d
ys = flsdata[:,j,k]

plt.plot(xs, ys ,'-o',lw=1,ms=3)

xs = flxlin1d
ys = flxdata[:,j*nvisu//(npoly+1),k*nvisu//(npoly+1)]
#ys = (flxdata[i*nvisu//(npoly+1),:,k*nvisu//(npoly+1)]+flxdata[i*nvisu//(npoly+1)+1,:,k*nvisu//(npoly+1)+1])/2

plt.plot(xs,ys,'g-', lw=1, ms=3)
#for x,y in zip(xs.reshape(-1,nvisu), ys.reshape(-1,nvisu)):
#    plt.plot(x,y,'g-', lw=1, ms=3)

plt.legend(['FLASH','FLEXI: nvisu = %d' % nvisu])
plt.savefig(outpath,bbox_inches='tight')
