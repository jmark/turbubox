#!/usr/bin/env pyturbubox

# stdlib
import os, sys, pickle
import numpy as np
from pathlib import Path

def put(msg):
    print(msg, file=sys.stderr, flush=True)

## ========================================================================= ##
## process commandline arguments

import argparse

pp = argparse.ArgumentParser(description = 'Batch Min Max Analysis Cache Files')

pp.add_argument(
    '--snkfp',
    help='path to store minmax.pickle',
    type=Path, required=True,
)

pp.add_argument(
    '--skip',
    type=bool,
    default=False,
    help='skip already produced files',
)

pp.add_argument(
    '--pickles',
    help='list of snapshot files',
    type=Path,nargs='*', required=True
)

ARGV = pp.parse_args()

key = 'pdf_vw'
subkey = 'dens'

dd = {
    key: {subkey: {'min': [999,999], 'max': [0,0]}},
}

for pid, srcfp in enumerate(ARGV.pickles):
    with open(str(srcfp), 'rb') as fh: data = pickle.load(fh)
    xs = data[key][subkey][1]
    ys = data[key][subkey][0]

    dd[key][subkey]['min'][1] = min(dd[key][subkey]['min'][1], np.min(xs))
    dd[key][subkey]['max'][1] = max(dd[key][subkey]['max'][1], np.max(xs))
                                                          
    dd[key][subkey]['min'][0] = min(dd[key][subkey]['min'][0], np.min(ys))
    dd[key][subkey]['max'][0] = max(dd[key][subkey]['max'][0], np.max(ys))

with open(str(ARGV.snkfp), 'wb') as fh: pickle.dump(dd, fh)
