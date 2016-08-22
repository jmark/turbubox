import itertools
import pandas as pd
import hashlib

def gen_key(unit,grid,solv,deli='/'):
    return '%s%s%d%s%s' % (unit,deli,grid,deli,solv)

def gen_amr_path(key):
    return '/srv/data/FLASH/stirturb/mach-2.0/%s' % key

def gen_pickle_path(key):
    cachedir = '/tmp'
    #uniqueid = hashlib.md5(key.encode()).hexdigest()
    uniqueid = key.replace('/','.')
    return '%s/time-evolution.%s.pandas.pickle' % (cachedir,uniqueid)

def gen_hdf5_path(key):
    cachedir = '/tmp'
    uniqueid = key.replace('/','.')
    return '%s/time-evolution.%s.pandas.h5' % (cachedir,uniqueid)

units = 'unit cgs'.split()
grids = [16,24,32,48,64]
solvs = '8w b3 b5 es'.split()

def load_datafiles(units,grids,solvs):
    keys = [gen_key(*key) for key in itertools.product(units,grids,solvs)]
    dfs  = [pd.read_hdf(gen_hdf5_path(key),'table') for key in keys]

    return dict(zip(keys,dfs))

def turntime(key):
    if 'unit' in key:
        return 1/2.0
    elif 'cgs' in key:
        return 4.385e+14
