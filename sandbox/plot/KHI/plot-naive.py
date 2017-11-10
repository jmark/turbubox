#!/usr/bin/env pyturbubox

import ulz
import sys
import cubicle
import argparse
import numpy as np

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 9})
import matplotlib.pyplot as plt

pp = argparse.ArgumentParser()

pp.add_argument('srcfile', help='path to couchdg snapshot file', type=Path)
pp.add_argument('snkfile', help='path to of png file', type=Path)
pp.add_argument('--title', type=ulz.URLhandler, default='')

ARGV = pp.parse_args()

fdata = cubicle.File(ARGV.srcfile)

dpi = 150

fig = plt.figure(figsize=(1920/dpi, 1080/dpi), dpi=dpi)

plt.xlim(fdata.domain.T[0,:])
plt.ylim(fdata.domain.T[1,:])

dens,velx,vely,velz,pres = fdata.get_prims(Nvisu=10,gamma=5/3)
carpet = dens[:,:,0].T
extent = fdata.domain.T.ravel()[0:4]

print(carpet.shape)

plt.imshow(
    carpet,
    extent = extent,
    vmin = 0.5, vmax = 1.2,
    origin='lower',
    interpolation = None,
    cmap = plt.get_cmap('cubehelix'),
)

cb = plt.colorbar(fraction=0.0456, pad=0.02, format='%1.2f')
cb.set_label('   pressure', labelpad=-40, x=1.15, y=1.045, rotation=0)

plt.title('{0}: t = {1: 4.2f}'.format(ARGV.title.strip(), fdata.time), y=1.01)

plt.savefig(
    str(ARGV.snkfile),
    bbox_inches='tight',
    dpi=dpi,
)
