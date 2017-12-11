#!/usr/bin/env pyturbubox

import box
import ulz
import pickle
import argparse
from pathlib import Path
import scipy.optimize

import numpy as np
from matplotlib import pylab as plt

pp = argparse.ArgumentParser(description = 'Plotting PDF')

pp.add_argument(
    '--destdir',
    help='path to store: <dir>/%%03d.png',
    type=Path, required=True,
)

pp.add_argument(
    '--profile',
    type=str,
    required=True,
)

pp.add_argument(
    '--skip',
    type=bool,
    default=False,
    help='skip already produced files',
)

pp.add_argument(
    '--parallel',
    help='enable parallel processes: 0 --> max. n procs, > 0 --> set n procs',
    type=int,
    default=-1,
)

pp.add_argument(
    '--pickles',
    help='list of pickle files',
    type=Path, nargs='*', required=True
)

ARGV = pp.parse_args()

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

def task(args):
    taskID, srcfp = args

    snkfp = ARGV.destdir / srcfp.with_suffix('.png').name
    if ARGV.skip and snkfp.exists() and snkfp.stat().st_mtime > srcfp.stat().st_mtime: return

    with open(str(srcfp), 'rb') as fh: data = box.Box(pickle.load(fh))
  
    if ARGV.profile == 'dens_vw':
        #xrange = yrange = None
        xrange = (-30. ,-22.0)
        yrange = (-4,0)
        key    = 'pdf_vw'
        subkey = 'dens'
        title  = 'Volume-weighted Density PDF at Sim.-Time t = {:1.2e}'.format(data.time)
        xlabel = r'log. scale density log$_{10}(\rho)$'
        ylabel = 'log. scale PDF'
        fitfun = gauss_density
        fitini = (-27,1)

    else:
        raise NotImplementedError("Unknown setup {}".format(ARGV.setup))

    fig = plt.figure(figsize=(12,6))

    if xrange: plt.xlim(*xrange)
    if yrange: plt.ylim(*yrange)

    xs = ulz.bins2xs(data[key][subkey][1])
    ys = data[key][subkey][0]

    xs = xs[ys>0]
    ys = ys[ys>0]

    ys,xs = ulz.moving_avg_1d(ys,xs,5)

    plt.plot(xs,np.log10(ys),lw=1.5, label='PDF')

    if True:
        try:
            coeff, var_matrix = scipy.optimize.curve_fit(fitfun, xs, ys, p0=fitini)
            _ys = fitfun(xs, *coeff)

            mu     = coeff[0]
            sigma  = coeff[1]
            mach   = sigma_mach_relation(2*sigma)

            print(str(srcfp), mu, sigma, mach, sep="\t")
            plt.plot(xs,np.log10(_ys),':', color='gray', label='gaussian fit: mu = {:.2f}, sigma = {:.4f}'.format(mu,sigma))
        except TypeError as e:
            print('Warning: ', str(e))

    plt.grid()
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()

    plt.tight_layout()
    plt.savefig(str(snkfp), format='png')

    print('Finished: ', str(snkfp))

# ========================================================================= ##
# process tasks in serial or parallel
if ARGV.parallel < 0:
    for i,x in enumerate(ARGV.pickles): task((i,x))
else:
    import multiprocessing as mpr
    nprocs = None if ARGV.parallel == 0 else ARGV.parallel
    mpr.Pool(nprocs,maxtasksperchild=1).map(task,enumerate(ARGV.pickles))
