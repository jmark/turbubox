#!/usr/bin/env pyturbubox

import box
import ulz
import pickle
import argparse
from pathlib import Path
import scipy.optimize

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

ARGV = pp.parse_args()

with open(str(ARGV.pickle), 'rb') as fh:
    data = box.Box(pickle.load(fh))
  
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

if ARGV.setup == 'dens_vw':
    #xrange = yrange = None
    xrange = (-29.5,-27.0)
    yrange = (-3,1)
    key    = 'pdf_vw'
    subkey = 'dens'
    title  = 'Volume-weighted Density PDF at Sim.-Time t = {:1.2e}'.format(data.time)
    xlabel = r'log. scale density log$_{10}(\rho)$'
    ylabel = 'log. scale PDF'
    fitfun = gauss_density
    fitini = (-27,1)

else:
    raise NotImplementedError("Setup '{}' is unknown.".format(ARGV.setup))

## ------------------------------------------------------------------------- #
fig = plt.figure(figsize=(12,6))

## ------------------------------------------------------------------------- #

if xrange:
    plt.xlim(*xrange)
if yrange:
    plt.ylim(*yrange)

xs = ulz.bins2xs(data[key][subkey][1])
ys = data[key][subkey][0]

xs = xs[ys>0]
ys = ys[ys>0]

ys,xs = ulz.moving_avg_1d(ys,xs,5)

plt.plot(xs,np.log10(ys),lw=1.5, label='PDF')

if True:
    coeff, var_matrix = scipy.optimize.curve_fit(fitfun, xs, ys, p0=fitini)
    _ys = fitfun(xs, *coeff)

    mu     = coeff[0]
    sigma  = coeff[1]
    mach   = sigma_mach_relation(2*sigma)

    print(str(ARGV.pickle), mu, sigma, mach, sep="\t")
    plt.plot(xs,np.log10(_ys),':', color='gray', label='gaussian fit: mu = {:.2f}, sigma = {:.4f}'.format(mu,sigma))

plt.grid()
plt.title(title)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.legend()

plt.tight_layout()
plt.savefig(str(ARGV.output), format='png')
