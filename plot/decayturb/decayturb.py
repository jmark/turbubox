#!/usr/bin/env pyturbubox

import pandas as pd
import numpy as np
import scipy.misc, scipy.ndimage, scipy.integrate
import flash, flexi, hopr, periodicbox, ulz, gausslobatto, interpolate #jmark
import glob,pickle
import functools
import matplotlib.ticker as mtick
import scipy.optimize
import datetime

from pathlib import Path

gamma = 5/3
root = Path('/home/jmark/projects/stirturb/data/final/decayturb/')
ttc  = 0.1

class Data:
    def __init__(self,root, ttc=1., label=None, color=None,line=None):
        self.root  = Path(root)
        self.ttc   = ttc
        self.label = label
        self.color = color
        self.line  = line
        
        self.load()
        
        self.postproc()
        
    def load(self):
        self.glob = self.load_globals()
        self.cool = self.load_cooling()
        self.forc = self.load_forcing()
        self.prof = self.load_profiling()
        self.anal = self.load_analysis()
    
    def load_profiling(self):
        fp = self.root / 'profiling.dat'
        names = 'time dt year month day tz hour minute sec msec'.split()
        df = pd.read_csv(fp, sep='\s+',names=names)
        df['datetime'] = df.apply(lambda row: datetime.datetime(*(int(x) for x in [
                                                row.year,row.month,row.day,row.hour,row.minute,row.sec,row.msec
                                            ])).timestamp(), axis=1)
        df['runtime'] = df.datetime - df.datetime[0]
        df['dyntime'] = df.time / self.ttc
        return df
    
    def load_analysis(self):
        fp = self.root / 'pickle'
        
        def load_pickle(fp):
            with open(fp, 'rb') as fh:
                return pickle.load(fh)

        pckls = [load_pickle(fp) for fp in sorted(glob.glob(str(fp) + '/*.pickle'))]
        
        for p in pckls:
            p['scalars'] = dict(time = p['time'], dyntime = p['time'] / ttc, taskid = p['taskid'])
            del p['time']
            del p['taskid']
        
        #keys    =  pckls[0].keys()
        keys = 'scalars mean msqu min max pws pws2 pdf'.split()
        subkeys = {key: pckls[0][key].keys() for key in keys}
        
        return {key: {subkey: np.array([p[key][subkey] for p in pckls]) for subkey in subkeys[key]} for key in keys}
    
    def postproc(self):
        pass

class DataFlexi(Data):
    def load_cooling(self):
        fp = self.root / 'cooling.dat'
        names = 'time dt ener_bf eint_bf ekin_bf ener_af eint_af ekin_af'.split()
        df = pd.read_csv(fp,sep='\s+',names=names)
        df['eint']  = df.eint_af
        df['dEint'] = df.ener_af - df.ener_bf
        df['dEint_dt'] = df.dEint / df.dt
        df['dyntime'] = df.time / self.ttc
        return df
    
    def load_forcing(self):
        fp = self.root / 'stirring.dat'
        names = 'time dt fv_dg ener eint ekin rmsv mach mach_max'.split()
        df = pd.read_csv(fp,sep='\s+',names=names)
        df['dyntime'] = df.time / self.ttc
        return df    

    def load_globals(self):
        return self.load_forcing()

class DataFlash(Data):
    def load_cooling(self):
        fp = self.root / 'polytrope.dat'
        names = 'time dt dEint dEint_dt'.split()
        df = pd.read_csv(fp,sep='\s+',names=names,skiprows=1)
        df['dyntime'] = df.time / self.ttc
        return df
    
    def load_globals(self):
        fp = self.root / 'flash.dat'
        names = 'time mass momx momy momz ener ekin eint emag'.split()
        df = pd.read_csv(fp,sep='\s+',names=names,skiprows=1)
        df['dyntime'] = df.time / ttc
        return df

    def load_forcing(self):
        fp = self.root / 'global_analysis.dat'
        names = 'time dt ekin dEkin dEkin_dt rmsv_vw rmsv'.split()
        df = pd.read_csv(fp,sep='\s+',names=names,skiprows=2)
        df['dyntime'] = df.time / self.ttc
        return df

if __name__ == '__main__':
    runs = dict( 
        eu_fv = DataFlexi(root / 'euler-fv', ttc=ttc,label='Euler Fin-Vo',color='blue',line='-'),
        eu_hy = DataFlexi(root / 'euler-hybrid',ttc=ttc,label='Euler Hybrid',color='blue',line='--'),

        mp_fv = DataFlexi(root / 'midpoint-fv',ttc=ttc,label='Midp. Fin-Vo',color='green',line='-'),
        mp_hy = DataFlexi(root / 'midpoint-hybrid',ttc=ttc,label='Midp. Hybrid',color='green',line='--'),

        rk_fv = DataFlexi(root / 'rk3-fv',ttc=ttc, label='RK-3  Fin-Vo',color='red',line='-'),
        rk_hy = DataFlexi(root / 'rk3-hybrid',ttc=ttc, label='RK-3  Hybrid',color='red',line='--'),
        
        bouc3 = DataFlash(root / 'bouchut3', ttc=ttc,label='Bouchut3',color='orange',line='-'),
        #bouc5 = DataFlash(root / 'bouchut5', ttc=ttc,label='Bouchut5'),
        bouc5 = None,
    )

    import pickle
    import sys

    with open(sys.argv[1], 'wb') as fh:
        pickle.dump(runs, fh)
