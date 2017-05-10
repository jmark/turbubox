#!/usr/bin/env pyturbubox

import pickle
import argparse
from pathlib import Path
from collections import namedtuple
from decayturb import *

pp = argparse.ArgumentParser(description = 'Plotting Powerspectra')

pp.add_argument(
    '--pickle',
    help='path of pickle file',
    type=Path, required=True,
)

pp.add_argument(
    '--png',
    help='path of png file',
    type=Path, required=True,
)

pp.add_argument(
    '--index',
    help='set index',
    type=int,
    required=True,
)

ARGV = pp.parse_args()

Runs = namedtuple('runs', 'eu_fv eu_hy mp_fv mp_hy rk_fv rk_hy bouc3 bouc5'.split())

with open(ARGV.pickle, 'rb') as fh:
    runs = Runs(**pickle.load(fh))

import numpy as np
from matplotlib import pylab as plt
from mpl_toolkits.mplot3d import axes3d
plt.rcParams['figure.figsize'] = (8, 8)

fig = plt.figure(figsize=(12,6))
axB = fig.add_subplot(111)
ii = ARGV.index
method = 'pws2'
key = 'ekin'

plt.title('Powerspectra: Total Kinetic Energy over Spatial Wave Number at Dynamic Time: t_d = %.1f' % (runs.eu_hy.anal['scalars']['time'][ii]/0.2),y=1.1)
#plt.ylim(12,20)

def plot(run):
    df = run.anal
    xs,ys = df[method][key][ii]
    #xs,ys = ulz.moving_avg_1d(xs.values,ys.values,7)
    xs = np.log10(xs)
    ys = np.log10(ys)
    axB.plot(xs,ys,label=run.label,lw=2,ls=run.line,color=run.color)
    
plot(runs.rk_fv)
plot(runs.rk_hy)

plot(runs.eu_fv)
plot(runs.eu_hy)

plot(runs.mp_fv)
plot(runs.mp_hy)

#plot(runs.bouc5)
plot(runs.bouc3)

if False:
    plt.axvline(x=a,color='black')
    plt.axvline(x=b,color='black')

    slope = np.mean(slopes)
    ofs = 18
    xs = np.linspace(a-0.1,b+0.1,10)
    plt.plot(xs,slope*xs + ofs,color='black',lw=2)

    plt.xticks(np.arange(0,2.8,0.2))

    plt.text(b+0.1,17.2,'<slopes> = %.2f' % slope, va='center', size='x-large')

xticks = np.arange(-0.6,3.0,0.2)

## top x-axis
axT = axB.twiny()
axT.set_xticks(xticks)
axT.set_xticklabels(['%.2f' % (10**x) for x in xticks])
axT.set_xlabel('spatial wave number')
axT.set_xlim(xticks[0],xticks[-1])

## bottom x-axis
axB.set_xlabel('log. scale spatial wave number')
axB.set_xlim(xticks[0],xticks[-1])
axB.set_xticks(xticks)
axB.grid()
axB.legend(ncol=1)

plt.tight_layout()
plt.savefig(str(ARGV.png), format='png')
