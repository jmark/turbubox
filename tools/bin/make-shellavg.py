#!/usr/bin/env pyturbubox

# stdlib
import os, sys, pickle
import numpy as np
import pathlib as pl

# jmark
import ulz
import couchdg
import shellavg

def PositiveInt(arg):
    x = int(arg)
    if x >= 0:
        return x
    else:
        raise ValueError("'%d' must be positive!" % x)

def log(msg):
    print(msg, file=sys.stderr)

## ========================================================================= ##
## process commandline arguments

import argparse

pp = argparse.ArgumentParser(description = 'Batch Produce Shellaverages')

pp.add_argument(
    '--nsamples',
    help='number samples taken by shellavg3d',
    type=PositiveInt,
)

pp.add_argument(
    'snapshot',
    help='list of snapshot files',
    type=pl.Path,
)

ARGV = pp.parse_args()

fdata = couchdg.Ribbon(ARGV.snapshot)

velx = fdata.stitch(1,7)
vely = fdata.stitch(2,7)

radii,avgs = shellavg.shell_avg_2d(np.sqrt(velx**2+vely**2), ARGV.nsamples)

print(np.sqrt(fdata.domsize[0]**2 + fdata.domsize[1]**2))

radii *= 0.5*np.sqrt(fdata.domsize[0]**2 + fdata.domsize[1]**2)/len(radii)

for radius, avg in zip(radii,avgs):
    print(radius, avg)
