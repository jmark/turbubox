#!/usr/bin/env pyturbubox

import pickle
import argparse
from pathlib import Path
from collections import namedtuple
from decayturb import *
import box
import re

import numpy as np
from matplotlib import pylab as plt

## ------------------------------------------------------------------------- #
pp = argparse.ArgumentParser(description = 'Plotting Powerspectrum')

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

def parseSlice(argstr):
    mm = re.match('(\d+):(\d+)', argstr)
    if mm:
        a,b = mm.groups()
        return slice(int(a),int(b))
    else:
        a = int(argstr)
        return slice(a,a+1)

pp.add_argument(
    '--index',
    help='set index',
    type=parseSlice, required=True,
)

pp.add_argument(
    '--key',
    help='keyword method: pws pws2 pws...',
    type=str, required=True,
)

pp.add_argument(
    '--subkey',
    help='keyword: dens rmsv ekin',
    type=str, required=True,
)

pp.add_argument('--ylabel')
pp.add_argument('--title')

pp.add_argument(
    '--xrange',
    type=lambda arg: tuple(float(x) for x in arg.split(',')),
    #default=(-4,4),
)

pp.add_argument(
    '--yrange',
    type=lambda arg: tuple(float(x) for x in arg.split(',')),
    #default=(-7,2),
)

pp.add_argument(
    '--fit'
)

ARGV = pp.parse_args()

with open(ARGV.pickle, 'rb') as fh:
    runs = box.Box(pickle.load(fh))
  
runs.order = [
    runs.bouc3, runs.bouc5, runs.flppm,
    runs.mp_fv, runs.mp_hy,
    runs.rk_fv, runs.rk_hy,
]

runs.bouc3.color = 'orange'
runs.bouc5.color = 'orange'
runs.flppm.color = 'purple'

runs.mp_fv.color = 'blue'
runs.mp_hy.color = 'blue'

runs.rk_fv.color = 'green'
runs.rk_hy.color = 'green'

## ------------------------------------------------------------------------- #
## define fitting functions

def task(run,ii):
    df = run.anal

    _xs,_ys, *slurp = df[ARGV.key][ARGV.subkey][ii][0]

    area = np.trapz(_ys,_xs)

    xs = np.linspace(np.log10(_xs[0]),np.log10(_xs[-1]),1024)
    
    # take average and deviation
    nys,cys = 0.,0.
    for _xs,_ys, *slurp in df[ARGV.key][ARGV.subkey][ii]:
        nys += 1.
        cys += np.interp(xs,np.log10(_xs),np.log10(_ys))
    ys = cys / nys

    # take average and deviation
    nys,cys = 0.,0.
    for _xs,_ys, *slurp in df[ARGV.key][ARGV.subkey][ii]:
        nys += 1.
        cys += (ys - np.interp(xs,np.log10(_xs),np.log10(_ys)))**2
    dys = np.sqrt(cys / nys)

    Xs.append(xs); Ys.append(ys)
    #ys,xs = ulz.moving_avg_1d(ys,xs,7)

    __xs = np.linspace(0.6,1.4,100)
    __ys = np.interp(__xs,xs,ys)
    slope, ofs = np.polyfit(__xs, __ys,1)

    if slurp:
        print(run.label, "\t", *['%.4f' % x for x in [slope, ofs, *slurp[0]]])
    else:
        print(run.label, "\t", *['%.4f' % x for x in [slope, ofs, area/512**6]])
 
    axB.fill_between(xs, ys-dys, ys+dys, facecolor='grey', alpha=0.5)
    axB.plot(xs,ys,label=run.label,lw=1.5,ls=run.line,color=run.color)

## ------------------------------------------------------------------------- #
## Setup figure

fig = plt.figure(figsize=(12,6))
axB = fig.add_subplot(111)

if ARGV.yrange is not None:
    plt.ylim(*ARGV.yrange)

if ARGV.ylabel:
    plt.ylabel(ARGV.ylabel)
else:
    plt.ylabel('log. scale FFT[f(k)]')

if ARGV.title:
    plt.title(ARGV.title, y=1.1)
else:
    plt.title('Powerspectra', y=1.1)
    #plt.title('Powerspectra: Total Kinetic Energy over Spatial Wave Number at Dynamic Time: t_d = %.1f' % (runs.eu_fv.anal['scalars']['dyntime'][ARGV.index][0]),y=1.1)

## ------------------------------------------------------------------------- #
Xs,Ys = list(),list()
for run in runs.order:
    task(run,ARGV.index)

## ------------------------------------------------------------------------- #
## Plot average of all curves
if False:
    Xs = np.mean(np.array(Xs),axis=0)
    Ys = np.mean(np.array(Ys),axis=0)
    plt.plot(Xs,Ys,':',label='average', color='black')

## ------------------------------------------------------------------------- #
if ARGV.fit is not None:
    xlim = plt.gca().get_xlim()
    ylim = plt.gca().get_ylim()

    f = gauss
    coeff, var_matrix = scipy.optimize.curve_fit(f, Xs, Ys, p0=[1., 0., 1.])
    xs = np.linspace(*xlim,100)
    ys = f(xs, *coeff)
    #label = 'gaussian fit: A=%6.3f, mu=%6.3f, sigma=%6.3f' % tuple(coeff)
    label = 'gaussian fit'
    plt.plot(xs,np.log10(ys),'--',label=label, color='black')

    f = skew_pdf
    #coeff = [2,1,0.7,-20]
    coeff = [*coeff,1.0]
    coeff, var_matrix = scipy.optimize.curve_fit(f, Xs, Ys, p0=coeff)
    xs = np.linspace(*xlim,100)
    ys = f(xs, *coeff)
    #label = 'pdf-skew fit: A=%6.3f, mu=%6.3f, sigma=%6.3f, alpha=%6.3f' % tuple(coeff)
    label = 'pdf-skew'
    plt.plot(xs,np.log10(ys),'-',label=label, color='black')

xticks = np.arange(-0.6,3.0,0.2)

_xticks = np.linspace(0.6,1.4,10)
axB.plot(_xticks, -1.1 * _xticks + 17.5, ls=':', lw=2, color='black', label='-1.1 • log10(k) + 17.5')
#axB.plot(_xticks, -19/9 * _xticks + 17.8, ls=':', lw=2, color='black', label='-19/9 • log10(k) + 17.8')
#axB.plot(_xticks, -5/3 * _xticks + 17, ls=':', lw=2, color='black', label='-5/3 • log10(k) + 17.8')

## top x-axis
axT = axB.twiny()
axT.set_xticks(xticks)
axT.set_xticklabels(['%.2f' % (10**x) for x in xticks])
axT.set_xlabel('spatial wave number k')
axT.set_xlim(xticks[0],xticks[-1])

## bottom x-axis
axB.set_xlabel('log. scale spatial wave number k')
axB.set_xlim(xticks[0],xticks[-1])
axB.set_xticks(xticks)
axB.grid()
axB.legend(ncol=1)

plt.tight_layout()
plt.savefig(str(ARGV.png), format='png')
