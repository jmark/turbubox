#!/usr/bin/env pyturbubox

# stdlib
import numpy as np
from pathlib import Path
import argparse

# jmark
import ulz
import cubicle
import shellavg

def PositiveInt(arg):
    x = int(arg)
    if x >= 0: return x
    else: raise ValueError("'%d' must be positive!" % x)

## ========================================================================= ##
## process commandline arguments

pp = argparse.ArgumentParser()

pp.add_argument(
    '--nvar',
    help='index number of variable, e.g.: 0 => dens, 1 => velx, ...',
    type=PositiveInt, required=True
)

pp.add_argument(
    '--nsamples',
    help='number samples taken by shellavg3d',
    type=PositiveInt,
)

pp.add_argument(
    '--meshfile',
    help='hopr mesh file for flexi snapshot files',
    type=Path,
)

pp.add_argument(
    '--gamma',
    help='number samples taken by shellavg3d',
    type=float,
)

pp.add_argument(
    'snapshot',
    help='snapshot file',
    type=Path,
)

ARGV = pp.parse_args()
cube = cubicle.File(ARGV.snapshot, meshfile=ARGV.meshfile)

if ARGV.nvar <= 3:
    box  = cube.as_box(ARGV.nvar)
else:
    prim = cube.get_prims(gamma=ARGV.gamma)
    box = prim[ARGV.nvar]

radii, avgs = shellavg.shell_avg(box, ARGV.nsamples)
radii *= np.sqrt(np.sum(cube.domsize**2))/np.sqrt(np.sum(np.array(box.shape)**2))

print('#', 'radius', 'shell-average')
for radius, avg in zip(radii,avgs):
    print(radius, avg)
