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

pp.add_argument('--srcfile', type=Path, required=True)
pp.add_argument('--snkfile', type=Path, required=True)
pp.add_argument('--profile', type=str.lower)
pp.add_argument('--title',   type=lambda x: ulz.URLhandler(x).strip())
pp.add_argument('--nvisu',   type=int)
pp.add_argument('--gamma',   type=float)

ARGV = pp.parse_args()

fdata = cubicle.File(ARGV.srcfile)
gamma = fdata.gamma if hasattr(fdata,'gamma') else ARGV.gamma

cbrange = (None,None)

dpi = 150
fig = plt.figure(figsize=(1920/dpi, 1080/dpi), dpi=dpi)

if fdata.domain.shape[1] > 2:
    dens,velx,vely,velz,pres = [x[:,:,0].T for x in fdata.get_prims(gamma=gamma)]
    extent = fdata.domain.T.ravel()[0:4]
else:
    dens,velx,vely,velz,pres = fdata.get_prims(Nvisu=ARGV.nvisu,gamma=gamma)
    extent = fdata.domain.T.ravel()

if ARGV.profile == 'density':
    carpet  = fdata.stitch(0,Nvisu=ARGV.nvisu)
    #cbrange = (0.1,30)
    #cbrange = (0.0,15)
    cblabel = '   density'
    cmap    = plt.get_cmap('cubehelix')

elif ARGV.profile == 'log10-density':
    carpet  = np.log10(dens)
    cbrange = (None,None)
    cbrange = (-2,2)
    cbrange = (-2,1.0)
    cblabel = ' log10 dens.'
    cmap    = plt.get_cmap('cubehelix')

elif ARGV.profile == 'pressure':
    carpet  = pres
    cbrange = (0.4,1.3)
    cbrange = (0.4,1.8)
    cblabel = '   pressure'
    cmap    = plt.get_cmap('cubehelix')

elif ARGV.profile == 'mach':
    carpet  = np.sqrt(velx**2+vely**2+velz**2)/np.sqrt(gamma*np.abs(pres/dens))
    cbrange = (0.0,1.0)
    cbrange = (0.0,1.8)
    cblabel = '    Mach'
    cmap    = plt.get_cmap('cubehelix')

elif ARGV.profile == 'vorticity':
    import vectoranalysis2D as v2
    carpet  = v2.curl(velx,vely,Dx=fdata.domsize[0]/velx.shape[0],Dy=fdata.domsize[1]/vely.shape[1])
    #cbrange = (-45,45)
    cbrange = (-50,50)
    #cbrange = (-15,15)
    cblabel = '  vorticity'
    cmap    = plt.get_cmap('gnuplot2')

elif ARGV.profile == 'moe':
    moe     = fdata.stitch(fdata.profiles['moe'],dname='profile')
    carpet  = np.log10(moe+1e-16)
    cbrange = (0.0,-5.0)
    cblabel = '   log10(moe)'
    cmap    = plt.get_cmap('gist_heat')

elif ARGV.profile == 'edof':
    carpet  = fdata.stitch(0,dname='edof')
    carpet  = np.log10(carpet-1e-16)
    cbrange = (-15.0,0.0)
    cblabel = '    edof'
    cmap    = plt.get_cmap('gist_heat')

elif ARGV.profile == 'blend':
    carpet  = fdata.stitch(0,dname='blend')
    carpet  = np.log10(carpet+1e-16)
    cbrange = (-5.0,0.0)
    cblabel = 'log10(blend)'
    cmap    = plt.get_cmap('gist_heat')

elif ARGV.profile == 'vpot':
    carpet = fdata.stitch(fdata.profiles['vpot'],dname='profile',Nvisu=ARGV.nvisu)
    cblabel = '   far-field pot.'
    cbrange = (-5,-25)
    cmap    = plt.get_cmap('gist_heat')

else:
    raise NotImplementedError('Unknown profile: ' + ARGV.profile)

plt.imshow(
    carpet.T,
    extent = np.roll(extent,2),
    vmin = cbrange[0], vmax = cbrange[1],
    origin='lower',
    interpolation = None,
    cmap = cmap,
)

cb = plt.colorbar(fraction=0.0456, pad=0.02, format='%1.2f')
cb.set_label(cblabel, labelpad=-40, x=1.15, y=1.045, rotation=0)

if ARGV.title:
    title = '{}: t = {: 6.3f}'.format(ARGV.title, fdata.time)
else:
    title = 't = {: 6.3f}'.format(fdata.time)

plt.title(title, y=1.01)
plt.savefig(str(ARGV.snkfile), bbox_inches='tight', dpi=dpi)
