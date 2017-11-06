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

fdens = fftshift(np.abs(rfftn(dens)))**2
fvels = fftshift(np.abs(rfftn(vels)))**2
fpres = fftshift(np.abs(rfftn(pres)))**2
fekin = fftshift(np.abs(rfftn(ekin)))**2

radii, shavg_dens = shell_avg(fdens, ARGV.nsamples, want_powerspectrum=True)
radii, shavg_vels = shell_avg(fvels, ARGV.nsamples, want_powerspectrum=True)
radii, shavg_pres = shell_avg(fpres, ARGV.nsamples, want_powerspectrum=True)
radii, shavg_ekin = shell_avg(fekin, ARGV.nsamples, want_powerspectrum=True)

print('#', 'radius', 'dens', 'vels', 'pres', 'ekin')
for radius, dens, vels, pres, ekin in zip(radii,shavg_dens,shavg_vels,shavg_pres,shavg_ekin):
    print(radius, dens, vels, pres, ekin)
