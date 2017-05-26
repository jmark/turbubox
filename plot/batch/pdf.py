#!/usr/bin/env pyturbubox

import pickle
import argparse
from pathlib import Path
from collections import namedtuple
from decayturb import *
import box
import re
import uncertainties as uc

import numpy as np
from matplotlib import pylab as plt

## ------------------------------------------------------------------------- #
pp = argparse.ArgumentParser(description = 'Plotting PDF')

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
    help='keyword method: pdf_vw pdf_mw',
    type=str, #required=True,
)

pp.add_argument(
    '--subkey',
    help='keyword: dens rmsv ekin',
    type=str, #required=True,
)

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
    '--fit', action='store_true',
)

ARGV = pp.parse_args()

#Runs = namedtuple('runs', 'eu_fv eu_hy mp_fv mp_hy rk_fv rk_hy bouc3 bouc5'.split())

with open(ARGV.pickle, 'rb') as fh:
    runs = box.Box(pickle.load(fh))
  
runs.order = [
    runs.bouc3, runs.bouc5, runs.flppm,
    runs.mp_fv, runs.mp_hy, runs.rk_fv, runs.rk_hy,
]

# colour table in HTML hex format
hexcols = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', 
           '#CC6677', '#882255', '#AA4499', '#661100', '#6699CC', '#AA4466',
           '#4477AA']

runs.bouc3.color = 'orange'
runs.bouc5.color = 'orange'
runs.flppm.color = 'purple'
runs.flppm.label = 'PPM     '

runs.mp_fv.color = 'blue'
runs.mp_hy.color = 'blue'

runs.rk_fv.color = 'green'
runs.rk_hy.color = 'green'

## ------------------------------------------------------------------------- #
## define fitting functions

def pdf(x):
    return np.exp(-x**2/2)

def cdf(x):
    return 0.5*(1+scipy.special.erf(x/np.sqrt(2)))

def skew_pdf(x,*p):
    A, mu, sigma, alpha = p
    x = (x-mu)/sigma
    return 2*A/sigma*pdf(x)*cdf(alpha*x)

def gauss(x,*p):
    A, mu, sigma = p
    x = (x-mu)/sigma
    return A*pdf(x)

def gauss_density(x,*p):
    mu, sigma = p
    x = (x-mu)/sigma
    return np.exp(-x**2/2)/np.sqrt(2*np.pi*sigma**2)

def sigma_mach_relation(sigma):
    zeta = 0.5 # ratio compressible to solenoidal forcing
    b = 1-2/3*zeta
    return np.sqrt(np.exp(sigma**2)-1)/b

def sigma_mach_deriv(sigma):
    zeta = 0.5 # ratio compressible to solenoidal forcing
    b = 1-2/3*zeta
    return sigma * np.exp(sigma**2) / np.sqrt(np.exp(sigma**2)-1) / b

## ------------------------------------------------------------------------- #
## Set custom configurations 

if ARGV.setup == 'default':
    index  = ARGV.index
    key    = ARGV.key
    subkey = ARGV.subkey
    xrange = ARGV.xrange
    yrange = ARGV.yrange
    title  = '%s/%s PDF at Dynamic Time: t_d = %.1f' % (
            key, subkey, 
            runs.order[0].anal['scalars']['dyntime'][index][0])
    xlabel = 'log. scale %s' % ARGV.subkey
    ylabel = 'log. scale PDF'
    fitfun = gauss
    fitini = [1,0,1]

elif ARGV.setup == 'density vw':
    index  = slice(20,40)
    key    = 'pdf_vw'
    subkey = 'dens'
    xrange = ARGV.xrange
    yrange = ARGV.yrange
    title  = 'Mass-weighted Density PDF between 2 and 4 Crossing Times'
    xlabel = 'log. scale density'
    ylabel = 'log. scale PDF'
    fitfun = gauss_density
    fitini = [0,1]

elif ARGV.setup == 'density mw':
    index  = slice(20,40)
    key    = 'pdf_mw'
    subkey = 'dens'
    xrange = ARGV.xrange
    yrange = ARGV.yrange
    title  = 'Mass-weighted Density PDF between 2 and 4 Crossing Times'
    xlabel = 'log. scale density'
    ylabel = 'log. scale PDF'
    fitfun = gauss_density
    fitini = [0,1]

else:
    raise NotImplementedError("Setup '%s' is unknown." % ARGV.setup)

## ------------------------------------------------------------------------- #
fig = plt.figure(figsize=(12,6))


## ------------------------------------------------------------------------- #
if xrange is not None:
    plt.xlim(*xrange)
if yrange is not None:
    plt.ylim(*yrange)

## ------------------------------------------------------------------------- #
## Plot average of all curves

Xs,Ys = list(),list()

for run in runs.order:
    ii = index
    df = run.anal

    # take average
    nbins,npdf = 0.,0.
    cbins,cpdf = 0.,0.
    for pdf,bins in df[key][subkey][ii]:
        cbins += bins
        nbins += 1.
        cpdf  += pdf
        npdf  += 1.

    bins = cbins / nbins
    pdf  = cpdf  / npdf

    # take average
    nbins,npdf = 0.,0.
    cbins,cpdf = 0.,0.
    for _pdf,_bins in df[key][subkey][ii]:
        cpdf  += (pdf - _pdf)**2
        npdf  += 1.

    dpdf = np.sqrt(cpdf/npdf)

    xs  = ulz.bins2xs(bins)
    ys  = pdf
    dys = dpdf
    
    #n  = 12
    #ys,xs = ulz.moving_avg_1d(ys,xs,n)
    #dys = ulz.moving_avg_1d(dys,n=n)
    
    Xs.append(xs)
    Ys.append(ys)
    
    dys = dys[ys>0]
    xs  = xs[ys>0]
    ys  = ys[ys>0]

    plt.fill_between(xs, np.log10(np.where(ys-dys > 0, ys-dys, np.nan)), np.log10(ys+dys), facecolor='grey', alpha=0.5)

    plt.plot(xs,np.log10(ys),lw=1.5,label=run.label,ls=run.line,color=run.color)
    #plt.plot(xs,ys,lw=1.5,label=run.label,ls=run.line,color=run.color)

    if ARGV.fit is not None:
        coeff, var_matrix = scipy.optimize.curve_fit(fitfun, xs, ys, p0=fitini)
        _ys = fitfun(xs, *coeff)
        #label = 'gaussian fit: A=%6.3f, mu=%6.3f, sigma=%6.3f' % tuple(coeff)
        label = 'fit'
        plt.plot(xs,np.log10(_ys),':',label='_nolegend_', color='gray')

        coeffL, var_matrixL = scipy.optimize.curve_fit(fitfun, xs, ys-dys, p0=fitini)
        coeffU, var_matrixU = scipy.optimize.curve_fit(fitfun, xs, ys+dys, p0=fitini)

        mu     = coeff[0]
        dmu    = np.abs(coeffL[0] - coeffU[0])

        sigma  = coeff[1]
        dsigma = np.abs(coeffL[1] - coeffU[1])

        mach   = sigma_mach_relation(2*sigma)
        dmach  = sigma_mach_deriv(2*sigma)*2*dsigma

        #plt.plot(xs,_ys,':',label='_nolegend_', color='gray')
        #print(run.label, mu,dmu,sigma,dsigma, mach,dmach)

        print(run.label, 
            uc.ufloat(mu,dmu), uc.ufloat(sigma,dsigma), uc.ufloat(mach,dmach),
            sep="\t")

if False:
    Xs = np.mean(np.array(Xs),axis=0)
    Ys = np.mean(np.array(Ys),axis=0)

    Xs = Xs[Ys>0]
    Ys = Ys[Ys>0]

    #plt.plot(Xs,np.log10(Ys),':',label='average', color='black')

    f = gauss
    coeff, var_matrix = scipy.optimize.curve_fit(f, Xs, Ys, p0=[1., 0., 1.])

    xlim = plt.gca().get_xlim()
    ylim = plt.gca().get_ylim()

    _xs = np.linspace(*xlim,100)
    _ys = f(_xs, *coeff)
    #label = 'gaussian fit: A=%6.3f, mu=%6.3f, sigma=%6.3f' % tuple(coeff)
    label = 'gaussian fit'
    plt.plot(_xs,np.log10(_ys),'--',label=label, color='black',lw=1)

    print('gaussian fit: ', *coeff)

if False and ARGV.fit is not None:
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


## ------------------------------------------------------------------------- #

plt.grid()
plt.legend(
    #loc='upper right'
)

plt.title(title)
plt.xlabel(xlabel)
plt.ylabel(ylabel)

plt.tight_layout()

if ARGV.output in 'show':
    plt.show() 
elif ARGV.output in 'none':
    pass
else:
    plt.savefig(str(ARGV.output), format='png')
