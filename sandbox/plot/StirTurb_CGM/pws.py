#!/usr/bin/env pyturbubox

# Following code is shit.

import pickle
import argparse
from pathlib import Path
from collections import namedtuple
import box
import re
import scipy.stats

import numpy as np
from matplotlib import pylab as plt

## ------------------------------------------------------------------------- #
pp = argparse.ArgumentParser(description = 'Plotting Powerspectrum')

pp.add_argument(
    '--setup',
    help='task identifier',
    type=str.lower, #required=True,
    default='default',
)

pp.add_argument(
    '--pickle',
    help='path of pickle file',
    type=Path, required=True,
)

pp.add_argument(
    '--output',
    help='path of png file',
    type=str, required=True,
)

ARGV = pp.parse_args()

with open(str(ARGV.pickle), 'rb') as fh:
    data = box.Box(pickle.load(fh))
  
## ------------------------------------------------------------------------- #
## Set custom configurations 

xlabelBot = r'log. scale spatial wave number log$_{10}(k)$'
xlabelTop = r'spatial wave number $k$'
ylabel = 'log. scale FFT[f(k)]'

if ARGV.setup == 'volume-weighted/velocity':
    key    = 'pws3' 
    subkey = 'pws3d_vw' 
    ylabel = r'log. scale Fourier transformed velocity log$_{10}(\hat{u})$'
    title  = r'Shell-averaged Powerspectra of Three-dimensional Volume-weighted' \
            + r' Velocity Field at Time $t_d$ = {:1.2e}'.format(data.time)

    yrange = (23,32)
    #xticks = np.arange(-0.6,2.2,0.2)
    #xticks = np.arange(-0.6,2.6,0.2)
    xticks = np.arange(-0.6,2.8,0.2)

    irange = (0.6,1.4)

else:
    raise NotImplementedError("Setup '%s' is unknown." % ARGV.setup)

## ------------------------------------------------------------------------- #
## Setup figure

fig = plt.figure(figsize=(12,6))
axB = fig.add_subplot(111)

if yrange:
    plt.ylim(*yrange)

plt.ylabel(ylabel)
plt.title(title, y=1.1)

## --------------------------------------------------------------------- #

plt.axvline(x=irange[0],ls=':', color='black')
plt.axvline(x=irange[1],ls=':', color='black')

xs = data[key][subkey]['m1'][0]
ys = data[key][subkey]['m1'][1]

xs = np.log10(xs)
ys = np.log10(ys)

axB.plot(xs,ys, label='numerical data')

## --------------------------------------------------------------------- #

__xs = np.linspace(*irange,100)
__ys = np.interp(__xs,xs,ys)
slope, ofs, r_value, p_value, stderr = scipy.stats.linregress(__xs, __ys)

line_xs = np.linspace(*plt.gca().get_xlim(),10)
axB.plot(line_xs, slope*line_xs + ofs+0.5, ls=':', lw=1, color='green', label='linear fit: slope = {:.2f}'.format(slope))

## --------------------------------------------------------------------- #

## top x-axis
axT = axB.twiny()
axT.set_xticks(xticks)
axT.set_xticklabels(['%.2f' % (10**x) for x in xticks])
axT.set_xlabel(xlabelTop)
axT.set_xlim(xticks[0],xticks[-1])

## bottom x-axis
axB.set_xlabel(xlabelBot)
axB.set_xlim(xticks[0],xticks[-1])
axB.set_xticks(xticks)
axB.grid()
axB.legend()

plt.tight_layout()
plt.savefig(str(ARGV.output), format='png')
print(str(ARGV.output))
