#!/usr/bin/env pyturbubox

import ulz
import sys
import couchdg
import cubicle
import argparse
import numpy as np

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 8})
import matplotlib.pyplot as plt

pp = argparse.ArgumentParser()

pp.add_argument('srcfile', help='path to couchdg snapshot file', type=Path)
pp.add_argument('snkfile', help='path to of png file', type=Path)
pp.add_argument('--title', type=ulz.URLhandler, default='')
pp.add_argument('--profile', type=str.lower, default='density')

ARGV = pp.parse_args()

fdata = cubicle.File(ARGV.srcfile)

dpi = 150
fig = plt.figure(figsize=(1920/dpi, 1080/dpi), dpi=dpi)

plt.xlim(fdata.domain.T[0,:])
plt.ylim(fdata.domain.T[1,:])

gamma = 5/3
if fdata.domain.shape[1] > 2:
    dens,velx,vely,velz,pres = [x[:,:,0].T for x in fdata.get_prims(Nvisu=10,gamma=gamma)]
    extent = fdata.domain.T.ravel()[0:4]
else:
    dens,velx,vely,velz,pres = fdata.get_prims(Nvisu=10,gamma=gamma)
    extent = fdata.domain.T.ravel()

if ARGV.profile == 'density':
    carpet  = dens
    cbrange = (0.5,1.2)
    cblabel = '   density'
    cmap    = plt.get_cmap('cubehelix')

if ARGV.profile == 'log10-density':
    carpet  = np.log10(dens)
    cbrange = (None,None)
    cbrange = (-2,2)
    cblabel = ' log10 dens.'
    cmap    = plt.get_cmap('cubehelix')

elif ARGV.profile == 'pressure':
    carpet  = pres
    cbrange = (0.4,1.3)
    cblabel = '   pressure'
    cmap    = plt.get_cmap('cubehelix')

elif ARGV.profile == 'mach':
    carpet  = np.sqrt(velx**2+vely**2+velz**2)/np.sqrt(gamma*np.abs(pres/dens))
    cbrange = (0.0,1.0)
    cblabel = '    Mach'
    cmap    = plt.get_cmap('cubehelix')

elif ARGV.profile == 'vorticity':
    import vectoranalysis2D as v2
    carpet  = v2.curl(velx,vely,Dx=fdata.domsize[0]/velx.shape[0],Dy=fdata.domsize[1]/vely.shape[1])
    cbrange = (-45,45)
    cblabel = '  vorticity'
    cmap    = plt.get_cmap('gnuplot2')

elif ARGV.profile == 'limiter':
    theta_min = fdata.stitch(0,dname='limiter')
    theta_max = fdata.stitch(1,dname='limiter')
    carpet  = np.where(theta_min < theta_max, theta_min, theta_max)
    cbrange = (0.98,1.0)
    cblabel = '   limiter'
    cmap    = plt.get_cmap('gist_heat')

else:
    raise NotImplementedError('Unknown profile: ' + ARGV.profile)

plt.imshow(
    carpet,
    extent = extent,
    vmin = cbrange[0], vmax = cbrange[1],
    origin='lower',
    interpolation = None,
    cmap = cmap,
)

cb = plt.colorbar(fraction=0.0456, pad=0.02, format='%1.2f')
cb.set_label(cblabel, labelpad=-40, x=1.15, y=1.045, rotation=0)

plt.title(ARGV.title.strip() + ": t = {: 6.3f}".format(fdata.time), y=1.01)

plt.savefig(str(ARGV.snkfile), bbox_inches='tight', dpi=dpi)
