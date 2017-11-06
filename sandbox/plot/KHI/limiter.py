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
matplotlib.rcParams.update({'font.size': 9})
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

theta_min = fdata.stitch(0)
theta_max = fdata.stitch(1)

carpet  = np.where(theta_min < theta_max, theta_min, theta_max)


extent  = fdata.domain.T.ravel()
cbrange = (0.98,1.0)
cblabel = '   Moe-Limiter'
cmap    = plt.get_cmap('cubehelix')

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
