#!/usr/bin/env pyturbubox

# stdlib
import os, sys, pickle
from pathlib import Path
import argparse

# 3rd party
import numpy as np
from numpy.fft import rfftn, fftshift

# turbubox
import ulz
import cubicle
from shellavg import shell_avg

pp = argparse.ArgumentParser(description = 'Batch Produce Powerspectra')

pp.add_argument(
    '--nsamples',
    help='number samples taken by shellavg',
    type=ulz.PositiveInt,
)

pp.add_argument(
    '--meshfile',
    help='hopr mesh file for flexi snapshot files',
    type=Path,
)

pp.add_argument(
    'snapshot',
    help='snapshot file',
    type=Path
)

ARGV = pp.parse_args()

cube = cubicle.File(ARGV.snapshot, meshfile=ARGV.meshfile)

dens, velx, vely, velz, pres = cube.get_prims()
vels = np.sqrt(velx**2+vely**2+velz**2)
ekin = dens/2 * (velx**2+vely**2+velz**2)

fdens = fftshift(np.abs(rfftn(dens)))
fvels = fftshift(np.abs(rfftn(vels)))
fpres = fftshift(np.abs(rfftn(pres)))
fekin = fftshift(np.abs(rfftn(ekin)))

radii, shavg_dens = shell_avg(fdens**2, ARGV.nsamples, mult_with_rsquare=True)
radii, shavg_vels = shell_avg(fvels**2, ARGV.nsamples, mult_with_rsquare=True)
radii, shavg_pres = shell_avg(fpres**2, ARGV.nsamples, mult_with_rsquare=True)
radii, shavg_ekin = shell_avg(fekin**2, ARGV.nsamples, mult_with_rsquare=True)

shavg_dens *= 4*np.pi
shavg_vels *= 4*np.pi
shavg_pres *= 4*np.pi
shavg_ekin *= 4*np.pi

print('#', 'radius', 'dens', 'vels', 'pres', 'ekin')
for radius, dens, vels, pres, ekin in zip(radii,shavg_dens,shavg_vels,shavg_pres,shavg_ekin):
    print(radius, dens, vels, pres, ekin)
