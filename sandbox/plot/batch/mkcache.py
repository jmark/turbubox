#!/usr/bin/env pyturbubox

import pickle
import sys

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

class Data:
    def __init__(self, id, root, ttc=1., label=None, color=None,line=None):
        self.id    = id
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
        self.fits = self.load_fits()
    
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
            p['scalars'] = dict(time = p['time'], dyntime = p['time'] / self.ttc, taskid = p['taskid'])
            del p['time']
            del p['taskid']
        
        keys = 'scalars mean msqu min max pws pws2 pdf_vw pdf_mw'.split()
        subkeys = {key: pckls[0][key].keys() for key in keys}
        
        retval = {key: {subkey: np.array([p[key][subkey] for p in pckls]) for subkey in subkeys[key]} for key in keys}

        for key in 'pws1d_vw pws3d_vw pws3d_mw'.split():
            ll = []
            for p in pckls:
                dd = p['pws3']
                ll.append([*dd[key]['m0'], dd[key]['areas']])
            retval['pws3/%s/m0' % key] = dict(rmsv = ll)

        for key in 'pws1d_vw pws3d_vw pws3d_mw'.split():
            ll = []
            for p in pckls:
                dd = p['pws3']
                ll.append([*dd[key]['m1'], dd[key]['areas']])
            retval['pws3/%s/m1' % key] = dict(rmsv = ll)

        return retval

    def load_fits(self):
        fp = self.root /'..'/'tgv/'
        return dict(
            vw_dens = pd.read_csv(str(fp / ('pdf/vw/density/%s.tgv' % self.id)), sep='\s+', names='id dyntime mu dmu sigma dsigma mach dmach'.split()),
            vw_vels = pd.read_csv(str(fp / ('pws/vw/velocity/%s.tgv' % self.id)), sep='\s+', names='id dyntime slope dslope ofs dofs area darea'.split()),
            mw_vels = pd.read_csv(str(fp / ('pws/mw/velocity/%s.tgv' % self.id)), sep='\s+', names='id dyntime slope dslope ofs dofs area darea'.split()),
            vw_ekin = pd.read_csv(str(fp / ('pws/vw/ekin/%s.tgv' % self.id)), sep='\s+', names='id dyntime slope dslope ofs dofs area darea areaSmall dareaSmall'.split()),
        )
 
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
        fp = self.root / 'forcing_analyze.dat'
        if not fp.exists():
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
        df['dyntime'] = df.time / self.ttc
        return df

    def load_forcing(self):
        fp = self.root / 'global_forcing.dat'
        if fp.exists():
            names = '''
                time dt
                EkinLastStp   EkinBefFor   EkinAftFor  dEkin_diss   dEkin_diss_dt  dEkin_forc   dEkin_forc_dt
                Vrms_vw_BF   Vrms_vw_AF   Vrms_mw_BF   Vrms_mw_AF
                Eint        dEint     dEint_dt
                Emag        dEmag     dEmag_dt
                Epot        dEpot     dEpot_dt
                dErad       dErad_dt
                E_bulk_BF   dE_bulk_corr   dE_bulk_corr_dt
                E_bulk_AF   dE_bulk_AF     dE_bulk_dt rmsv_af rmsv_bf
            '''.split()

            df = pd.read_csv(fp,sep='\s+',names=names,skiprows=2)
            df['ekin'] = df.EkinBefFor
            df['rmsv'] = df.Vrms_mw_AF
            df['eint'] = df.Eint
            df['ener'] = df.eint + df.ekin
            df['dyntime'] = df.time / self.ttc
     
            return df

        fp = self.root / 'global_analysis.dat'
        if fp.exists():
            names = 'time dt Ekin dEkin dEkin_dt rmsv_vw rmsv'.split()
            df = pd.read_csv(fp,sep='\s+',names=names,skiprows=2)
            df['dyntime'] = df.time / self.ttc

            return df
        
        raise OSError('%s not found.' % str(fp))
        
if __name__ == '__main__':

    import argparse
    pp = argparse.ArgumentParser(description = 'Plotting PDFs')

    pp.add_argument(
        '--output',
        help='path where to save pickle cache file',
        type=Path, required=True,
    )

    pp.add_argument(
        '--root',
        help='root path where to find data',
        type=Path, required=True,
    )

    pp.add_argument(
        '--ttc',
        help='crossing time',
        type=float, required=True,
    )

    ARGV = pp.parse_args()

    runs = dict(
        eu_fv = DataFlexi('eu_fv', ARGV.root / 'euler-fv', ttc=ARGV.ttc,label='Euler Fin-Vo',color='cyan',line='-'),
        #eu_hy = DataFlexi('eu_hy', ARGV.root / 'euler-hybrid',ttc=ARGV.ttc,label='Euler Hybrid',color='blue',line='--'),

        mp_fv = DataFlexi('mp_fv', ARGV.root / 'midpoint-fv',ttc=ARGV.ttc,label='Midp. Fin-Vo',color='blue',line='-'),
        mp_hy = DataFlexi('mp_hy', ARGV.root / 'midpoint-hybrid',ttc=ARGV.ttc,label='Midp. Hybrid',color='blue',line='--'),

        rk_fv = DataFlexi('rk_fv', ARGV.root / 'rk3-fv',ttc=ARGV.ttc, label='RK-3  Fin-Vo',color='green',line='-'),
        rk_hy = DataFlexi('rk_hy', ARGV.root / 'rk3-hybrid',ttc=ARGV.ttc, label='RK-3  Hybrid',color='green',line='--'),

        bouc3 = DataFlash('bouc3', ARGV.root / 'bouchut3', ttc=ARGV.ttc,label='Bouchut3',color='orange',line='-'),
        bouc5 = DataFlash('bouc5', ARGV.root / 'bouchut5', ttc=ARGV.ttc,label='Bouchut5',color='orange',line='--'),

        flppm = DataFlash('flppm', ARGV.root / 'ppm', ttc=ARGV.ttc,label='PPM',color='purple',line='-'),
    )

    runs['bouc3'].prof.runtime *= 512/768
    runs['bouc5'].prof.runtime *= 512/768
    runs['flppm'].prof.runtime *= 512/768

    runs['order'] = [runs[x] for x in 'bouc3 bouc5 flppm eu_fv mp_fv mp_hy rk_fv rk_hy'.split()]

    with open(str(ARGV.output), 'wb') as fh:
        pickle.dump(runs, fh)
