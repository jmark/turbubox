#!/usr/bin/env python3

import numpy as np
import pandas as pd

import glob

import itertools
import os
import multiprocessing

from pyflash import FlashFile

from ulz import *

import numba

@numba.jit(nopython=True,target='cpu')
def powerspectrum(fekin):

    (Nx,Ny,Nz) = fekin.shape
    nws = np.zeros(int(np.sqrt(Nx**2+Ny**2+Nz**2)))
    pws = np.zeros(int(np.sqrt(Nx**2+Ny**2+Nz**2)))

    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                I = i - Nx/2
                J = j - Ny/2
                K = k - Nz/2
                r = np.sqrt(I**2+J**2+K**2)
                pws[int(r)] += r**2 * fekin[i,j,k]
                nws[int(r)] += 1

    for i,n in enumerate(nws):
        if n > 0:
            pws[i] /= n

    return pws                

def calcTimeEvolution(path):
    fpaths = sorted(glob.glob('%s/amr/flash_hdf5_chk_*' % path))
    reap = dict([(key,[]) for key in 'time step ekin power'.split()])

    for i,fp in enumerate(fpaths):
        flash = FlashFile(fp)

        # scalar
        time = flash.meta['real scalars']['time']
        step = flash.meta['integer scalars']['nstep']
        c_s  = flash.meta['real runtime']['c_ambient']
        rho0 = flash.meta['real runtime']['rho_ambient']
        Delta = flash.meta['spacing'][0]
        
        # 3d
        rho = flash.get_box('dens')
        vel = [flash.get_box('vel'+dim) for dim in 'x y z'.split()]

        vol  = rho.shape[0]*rho.shape[1]*rho.shape[2]
        ekin = rho/2 * (vel[0]**2 + vel[1]**2 + vel[2]**2) / rho0 / vol / c_s**2

        from numpy.fft import rfftn, fftshift
        fekin = fftshift(np.abs(rfftn(ekin)))**2

        pws = powerspectrum(fekin)

        # collect data we want
        reap['time'].append(time)
        reap['step'].append(step)
        reap['ekin'].append(ekin)
        reap['power'].append(pws)
        
        #if i % 20 == 0:
        #    print(i,end=" ",flush=True)

    return pd.DataFrame.from_dict(reap)

def gen_datafile(args):
    unit,grid,solv = args

    key     = gen_key(unit,grid,solv)
    srcpath = gen_amr_path(key) # source
    snkpath = gen_hdf5_path(key)  # sink

    df = calcTimeEvolution(srcpath)
    df.to_hdf(snkpath,'table_power',append=True) 

    print(key,flush=True)

units = 'cgs'.split()
grids = [64]
solvs = '8w b3'.split()

params = itertools.product(units,grids,solvs)

pool = multiprocessing.Pool(6)
stat = pool.map(gen_datafile, params)
pool.close()

# for param in params:
#     print(gen_key(*param),end=' => ',flush=True)
#     gen_datafile(param)    
#     print('')
