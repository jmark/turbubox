#!/usr/bin/env pyturbubox

# stdlib
import os, sys, pickle
import numpy as np
from numpy.fft import fftn, fftshift
import pathlib as pl
import scipy.interpolate

# jmark
import periodicbox, ulz
from shellavg import shell_avg_3d
from defer_signals import DeferSignals

def PositiveInt(arg):
    x = int(arg)
    if x >= 0:
        return x
    else:
        raise ValueError("'%d' must be positive!" % x)

def log(msg):
    print(msg, file=sys.stderr, flush=True)

## ========================================================================= ##
## process commandline arguments

import argparse

pp = argparse.ArgumentParser(description = 'Batch Produce Powerspectra')

pp.add_argument(
    '--destdir',
    help='path to store: <dir>/%%03d.pickle',
    type=pl.Path, required=True,
)

pp.add_argument(
    '--parallel',
    help='enable parallel processes: -1 --> serial, 0 --> max. procs, N > 0 --> set n procs',
    type=int,
    default=-1,
)

pp.add_argument(
    '--gamma',
    help='set gamma (default: 5/3)',
    type=float,
    default=5./3,
)

pp.add_argument(
    '--skip',
    help='skip already produced files',
    action='store_true',
)

pp.add_argument(
    '--nsamples',
    help='number samples taken by shellavg3d',
    type=PositiveInt,
)

pp.add_argument(
    'snapshots',
    help='list of snapshot files',
    type=pl.Path,nargs='*',
)

ARGV = pp.parse_args()

## ========================================================================= ##
## Setup environment

if not ARGV.destdir.exists():
    ARGV.destdir.mkdir()

## ========================================================================= ##
## Routines

def powerspectrum(data, nshells=1024, npoints=512):
    def body_centered_linspace(infimum, supremum, nNodes):
        domsize  = np.abs(supremum - infimum)
        offset   = domsize / nNodes / 2.
        return np.linspace(infimum  + offset, supremum - offset, nNodes, endpoint=True)

    nx,ny,nz = data.shape
    
    # prepare interpolator
    bclx = body_centered_linspace(-nx/2,nx/2,nx)
    bcly = body_centered_linspace(-ny/2,ny/2,ny)
    bclz = body_centered_linspace(-nz/2,nz/2,nz)
    interpolate = scipy.interpolate.RegularGridInterpolator(
        (bclx,bcly,bclz), data, bounds_error=False, fill_value=np.NaN, method='linear')

    # array of radiis
    r_max = np.sqrt(nx**2+ny**2+nz**2)
    radii = body_centered_linspace(0, r_max, nshells)
    
    # spherical coordinates
    radiis = np.repeat(radii, npoints).reshape(-1,npoints)
    radiis += (np.random.rand(*radiis.shape) - 0.5) * r_max/len(radii)
    thetas =   np.pi*np.random.rand(nshells,npoints)
    phis   = 2*np.pi*np.random.rand(nshells,npoints)

    # cartesian coordinates
    xs = radiis * np.sin(thetas) * np.cos(phis)
    ys = radiis * np.sin(thetas) * np.sin(phis)
    zs = radiis * np.cos(thetas)

    # array of random points on shells in cartesian coordinates
    shells = np.array([xs,ys,zs]).transpose(1,2,0).reshape(-1,3)
    
    # interpolate shell points, trim NaNs and average
    vals = interpolate(shells).reshape(-1,npoints)
    NaNs = np.isnan(vals); vals[NaNs] = 0
    weights = npoints - np.sum(NaNs,axis=1)
    indices = np.where(weights > 0)

    avgs = np.sum((radiis**2*vals)[indices],axis=1) / weights[indices]
    
    return radii[indices], 4*np.pi*avgs

def pws(data):
    return shell_avg_3d( fftshift(np.abs(fftn(data)))**2, ARGV.nsamples )

def pws2(data):
    return powerspectrum( fftshift(np.abs(fftn(data)))**2, nshells = 3*len(data) )

def pdf(data):
    try:
        return np.histogram(data, bins=1024*4, range=(np.nanmin(data),np.nanmax(data)), density=True)
    except ValueError:
        return np.histogram(data, bins=1024*4, range=(-1,1), density=True)

def pws1d_vw(dens,velx,vely,velz,pres):
    fvelx = np.fft.fftn(velx)
    #fvely = np.fft.fftn(vely)
    #fvelz = np.fft.fftn(velz)

    fuck = 0.5*np.abs(fvelx)**2 #+ np.abs(fvely)**2 + np.abs(fvelz)**2)
    fuck = np.fft.fftshift(fuck)

    A0 = 0.5*np.mean(velx**2) #+ 0.5*np.mean(vely**2) + 0.5*np.mean(velz**2))
    A1 = np.sum(fuck) / len(fuck)**6

    # method 01
    radii,shells = shell_avg_3d(fuck)
    A2 = np.trapz(shells, radii) / len(fuck)**6

    # method 02
    radii2,shells2 = powerspectrum(fuck,nshells=3*len(fuck))
    A3 = np.trapz(shells2, radii2) / len(fuck)**6

    return dict(areas = [A0,A1,A2,A3], m0 = [radii,shells], m1 = [radii2,shells2])


def pws3d_vw(dens,velx,vely,velz,pres):
    fvelx = np.fft.fftn(velx)
    fvely = np.fft.fftn(vely)
    fvelz = np.fft.fftn(velz)

    fuck = 0.5*np.abs(fvelx)**2 + 0.5*np.abs(fvely)**2 + 0.5*np.abs(fvelz)**2
    fuck = np.fft.fftshift(fuck)

    A0 = 0.5*np.mean(velx**2) + 0.5*np.mean(vely**2) + 0.5*np.mean(velz**2)
    A1 = np.sum(fuck) / len(fuck)**6

    # method 01
    radii,shells = shell_avg_3d(fuck)
    A2 = np.trapz(shells, radii) / len(fuck)**6

    # method 02
    radii2,shells2 = powerspectrum(fuck,nshells=3*len(fuck))
    A3 = np.trapz(shells2, radii2) / len(fuck)**6

    return dict(areas = [A0,A1,A2,A3], m0 = [radii,shells], m1 = [radii2,shells2])

def pws3d_mw(dens,velx,vely,velz,pres):
    fvelx = np.fft.fftn(dens**0.5*velx)
    fvely = np.fft.fftn(dens**0.5*vely)
    fvelz = np.fft.fftn(dens**0.5*velz)

    fuck = 0.5*(np.abs(fvelx)**2 + np.abs(fvely)**2 + np.abs(fvelz)**2)
    fuck = np.fft.fftshift(fuck)

    A0 = 0.5*(np.mean((dens**0.5*velx)**2) + np.mean((dens**0.5*vely)**2) + np.mean((dens**0.5*velz)**2))
    A1 = np.sum(fuck) / len(fuck)**6

    # method 01
    radii,shells = shell_avg_3d(fuck)
    A2 = np.trapz(shells, radii) / len(fuck)**6

    # method 02
    radii2,shells2 = powerspectrum(fuck,nshells=3*len(fuck))
    A3 = np.trapz(shells2, radii2) / len(fuck)**6

    return dict(areas = [A0,A1,A2,A3], m0 = [radii,shells], m1 = [radii2,shells2])

def analysis(taskid, srcfp):
    box = periodicbox.File(srcfp, mode='r')
    dens, velx, vely, velz, pres = box.get_prims()
    gamma = ARGV.gamma

    rmsv = np.sqrt(velx**2+vely**2+velz**2) # rmsv
    rhov = dens**(1/3) * rmsv
    ekin = dens/2. * rmsv**2
    mach = rmsv / np.sqrt(gamma*pres/dens)

    fv = box.domainsize / np.array(dens.shape) # finite volume dimensions
    vort = np.mean(fv)**5/12.0 * np.abs(dens) * ulz.norm(*ulz.curl(velx,vely,velz,fv[0],fv[1],fv[2]))

    keys  = 'dens  pres  rmsv  ekin  mach  vort  rhov'.split()
    vals  = [dens, pres, rmsv, ekin, mach, vort, rhov]
    apply = lambda f: {k: f(v) for k,v in zip(keys,vals)}

    return dict(
        taskid = taskid,
        time   = box.time,

        mean = apply(np.mean),
        msqu = apply(lambda dd: np.mean(dd**2)),

        min  = apply(np.min),
        max  = apply(np.max),

        pws  = apply(pws),
        pws2 = apply(pws2),

        pdf_vw  = apply(lambda x: pdf(np.log10(x))),
        pdf_mw  = apply(lambda x: pdf(np.log10(dens*x/np.sum(dens)))),

        pws3 = dict(
            pws1d_vw = pws1d_vw(dens,velx,vely,velz,pres),
            pws3d_vw = pws3d_vw(dens,velx,vely,velz,pres),
            pws3d_mw = pws3d_mw(dens,velx,vely,velz,pres),
        ),
    )

## ========================================================================= ##
## prepare task

def task(taskid, srcfp):
    snkfp = ARGV.destdir / srcfp.with_suffix('.pickle').name

    if ARGV.skip and snkfp.exists() and snkfp.stat().st_mtime > srcfp.stat().st_mtime:
        return
    else:
        retval = analysis(taskid, srcfp)

        # prevent inconsistent states
        tmpfp  = snkfp.with_suffix('.tmp')
        with open(str(tmpfp), 'wb') as fd:
            pickle.dump(retval, fd)
        tmpfp.replace(snkfp) 

        log(str(snkfp))

## ========================================================================= ##
## do work

if ARGV.parallel >= 0:
    import multiprocessing as mpr
    def _task(args):
        return task(*args)
    nprocs = None if ARGV.parallel == 0 else ARGV.parallel
    mpr.Pool(nprocs,maxtasksperchild=1).map(_task,enumerate(ARGV.snapshots))
else:
    for i,x in enumerate(ARGV.snapshots):
        task(i,x)
