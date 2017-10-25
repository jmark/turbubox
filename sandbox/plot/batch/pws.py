#!/usr/bin/env pyturbubox

# Following code is shit.

import pickle
import argparse
from pathlib import Path
from collections import namedtuple
from decayturb import *
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
    type=parseSlice, #required=True,
)

pp.add_argument(
    '--key',
    help='keyword method: pws pws2 pws...',
    type=str, #required=True,
)

pp.add_argument(
    '--subkey',
    help='keyword: dens rmsv ekin',
    type=str, #required=True,
)

pp.add_argument('--ylabel')
pp.add_argument('--title')

pp.add_argument(
    '--xrange',
    type=lambda arg: tuple(float(x) for x in arg.split(':')),
    #default=(-4,4),
)

pp.add_argument(
    '--yrange',
    type=lambda arg: tuple(float(x) for x in arg.split(':')),
    #default=(-7,2),
)

pp.add_argument(
    '--fit', action='store_true',
)

ARGV = pp.parse_args()

with open(ARGV.pickle, 'rb') as fh:
    runs = box.Box(pickle.load(fh))
  
## ------------------------------------------------------------------------- #
## Set custom configurations 

# setting default values
output = ARGV.output
index  = ARGV.index
key    = ARGV.key
subkey = ARGV.subkey
xrange = ARGV.xrange
yrange = ARGV.yrange
title  = '%s/%s powerspectra at Dynamic Time: t_d = %.2f' % (
        key, subkey, 
        runs.order[0].anal['scalars']['dyntime'][index][0])
xlabelBot = r'log. scale spatial wave number log$_{10}(k)$'
xlabelTop = r'spatial wave number $k$'
ylabel = 'log. scale FFT[f(k)]'
fit = ARGV.fit

dyntime = runs.order[0].anal['scalars']['dyntime'][index][0]

if False:
    pass

elif ARGV.setup == 'mass-weighted/velocity':
    key    = 'pws3/pws3d_mw/m1' 
    subkey = 'rmsv' 
    ylabel = r'log. scale Fourier transformed velocity log$_{10}(\hat{u})$'
    title  = r'Shell-averaged Powerspectra of Three-dimensional Mass-weighted' \
            + r' Velocity Field at Dynamic Time $t_d$ = %.1f' % (dyntime)

elif ARGV.setup == 'volume-weighted/velocity':
    key    = 'pws3/pws3d_vw/m1' 
    subkey = 'rmsv' 
    ylabel = r'log. scale Fourier transformed velocity log$_{10}(\hat{u})$'
    title  = r'Shell-averaged Powerspectra of Three-dimensional Volume-weighted' \
            + r' Velocity Field at Dynamic Time $t_d$ = %.1f' % (dyntime)

elif ARGV.setup == 'volume-weighted/density':
    key    = 'pws2' 
    subkey = 'dens' 
    ylabel = r'log. scale Fourier transformed density log$_{10}(\hat{\rho})$'
    title  = r'Shell-averaged Powerspectra of Volume-weighted' \
            + r' Density Field at Dynamic Time $t_d$ = %.1f' % (dyntime)

elif ARGV.setup == 'volume-weighted/ekin':
    key    = 'pws2' 
    subkey = 'ekin' 
    ylabel = r'log. scale Fourier transformed kinetic energy log$_{10}(\hat{\mathcal{K}})$'
    title  = r'Shell-averaged Powerspectra of Volume-weighted' \
            + r' Kinetic Energy Field at Dynamic Time $t_d$ = %.1f' % (dyntime)

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

## ------------------------------------------------------------------------- #
Xs,Ys = list(),list()

ii = index
for run in runs.order:
    df = run.anal[key][subkey][ii]

    # take averaged area
    nn, carea = 0.,0.
    for _xs,_ys, *slurp in df:
        nn += 1.
        carea += np.trapz(_ys,_xs) / 512**6
    area = carea / nn

    # take averaged area deviation
    nn, carea = 0.,0.
    for _xs,_ys, *slurp in df:
        nn += 1.
        carea += (area - np.trapz(_ys,_xs) / 512**6)**2
    darea = 0.1 * area + np.sqrt(carea / nn)

    ## --------------------------------------------------------------------- #

    # get small scale area
    # set analysis domain
    _xs_ = np.linspace(64,256,512)

    exactAs = []
    # take averaged area
    nn, carea = 0.,0.
    for _xs,_ys, *slurp in df:
        _ys_ = np.interp(_xs_,_xs,_ys)
        nn += 1.
        carea += np.trapz(_ys_,_xs_) / 512**6
        exactAs.append(slurp)
    areaSmall = carea / nn

    exactAs = np.array(exactAs).T

    # take averaged area deviation
    nn, carea = 0.,0.
    for _xs,_ys, *slurp in df:
        _ys_ = np.interp(_xs_,_xs,_ys)
        nn += 1.
        carea += (areaSmall - np.trapz(_ys_,_xs_) / 512**6)**2
    dareaSmall = 0.1 * areaSmall + np.sqrt(carea / nn)

    ## --------------------------------------------------------------------- #

    # get original data
    _xs,_ys, *slurp = df[0]

    # set analysis domain
    xs = np.linspace(np.log10(_xs[0]),np.log10(_xs[-1]),1024)
    
    # take averaged codomain
    nys,cys = 0.,0.
    for _xs,_ys, *slurp in df:
        nys += 1.
        cys += np.interp(xs,np.log10(_xs),np.log10(_ys))
    ys = cys / nys

    # take averaged deviation
    nys,cys = 0.,0.
    for _xs,_ys, *slurp in df:
        nys += 1.
        cys += (ys - np.interp(xs,np.log10(_xs),np.log10(_ys)))**2
    dys = 1.4 * np.sqrt(cys / nys)

    Xs.append(xs); Ys.append(ys)
    #ys,xs = ulz.moving_avg_1d(ys,xs,7)

    ## --------------------------------------------------------------------- #

    if fit:
        __xs = np.linspace(0.6,1.2,100)

        __ys = np.interp(__xs,xs,ys)
        #slope, ofs   = np.polyfit(__xs, __ys,1)
        slope, ofs, r_value, p_value, stderr = scipy.stats.linregress(__xs, __ys)

        __ys = np.interp(__xs,xs,ys-dys)
        #slopeL, ofsL = np.polyfit(__xs, __ys,1)
        slopeL, ofsL, r_valueL, p_valueL, stderrL = scipy.stats.linregress(__xs, __ys)

        __ys = np.interp(__xs,xs,ys+dys)
        slopeU, ofsU = np.polyfit(__xs, __ys,1)
        slopeU, ofsU, r_valueU, p_valueU, stderrU = scipy.stats.linregress(__xs, __ys)

        dslope = np.abs(slopeU - slopeL) + stderr + stderrL + stderrU
        dofs   = np.abs(ofsU - ofsL) + stderr + stderrL + stderrU 

        #print(run.label, "\t", *['%.4f' % x for x in [slope, ofs, area/512**6]])
        #print(run.label, "\t", *['%.4f' % x for x in [slope, ofs, dslope, dofs]])
        #print(run.id, dyntime, 
        #    uc.ufloat(slope, dslope),
        #    uc.ufloat(ofs, dofs),
        #    uc.ufloat(area, darea),
        #    sep="\t")

        #print(exactAs)
        #sys.exit()

        foo = sum([[np.mean(x[0]), np.std(x[0])] for x in exactAs],[])
        #print(run.id, dyntime, slope, dslope, ofs, dofs, area, darea, areaSmall, dareaSmall, sep="\t")
        print(run.id, *foo, (foo[4] + foo[6])/2., (foo[5] + foo[7])/2., sep="\t")

        line_xs = np.linspace(*plt.gca().get_xlim(),10)
        axB.plot(line_xs, slope*line_xs + ofs, ls=':', lw=1, color=run.color, label='_nolegend_')
 
    ## --------------------------------------------------------------------- #

    axB.fill_between(xs, ys-dys, ys+dys, facecolor='grey', alpha=0.5)
    axB.plot(xs,ys,label=run.label,lw=1.5,ls=run.line,color=run.color)

## ------------------------------------------------------------------------- #
## Plot average of all curves
if False:
    Xs = np.mean(np.array(Xs),axis=0)
    Ys = np.mean(np.array(Ys),axis=0)
    plt.plot(Xs,Ys,':',label='average', color='black')

## ------------------------------------------------------------------------- #
if False and fit:
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

if False and fit:
    _xticks = np.linspace(0.6,1.4,10)
    axB.plot(_xticks, -1.1 * _xticks + 17.5, ls=':', lw=2, color='black', label='-1.1 • log10(k) + 17.5')
    #axB.plot(_xticks, -19/9 * _xticks + 17.8, ls=':', lw=2, color='black', label='-19/9 • log10(k) + 17.8')
    #axB.plot(_xticks, -5/3 * _xticks + 17, ls=':', lw=2, color='black', label='-5/3 • log10(k) + 17.8')

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
axB.legend(ncol=1)

plt.tight_layout()

if output in 'show':
    plt.show() 
elif output in 'none':
    pass
else:
    plt.savefig(str(output), format='png')
