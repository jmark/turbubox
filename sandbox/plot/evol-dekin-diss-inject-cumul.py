#!/usr/bin/env python3

import sys
import pickle
import numpy as np
from matplotlib import pylab as plt
import matplotlib
matplotlib.rc('font', family='DejaVu Sans')
matplotlib.rcParams.update({'font.size': 20})

def moving_sum(xs, n=2):
    return np.array([np.sum(x) for x in xs[:xs.shape[0]//n * n].reshape(-1,n)])

def moving_avg(xs, n=2):
    return np.array([np.mean(x) for x in xs[:xs.shape[0]//n * n].reshape(-1,n)])

import cycler
ccycle = cycler.cycler('color', ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00'])
plt.rc('axes', prop_cycle=ccycle)

# [01]time          
# [02]dt  
# [03]EkinLastStp   
# [04]EkinBefFor    [05]EkinAftFor
# [06]dEkin(diss)   [07] ([06]/dt)  
# [08]dEkin(forc)   [09] ([08]/dt)
# [10]Vrms_vw_BF    [11]Vrms_vw_AF      [12]Vrms_mw_BF  [13]Vrms_mw_AF
# [14]Eint          [15]dEint           [16]dEint/dt    [17]Emag            [18]dEmag
# [19]dEmag/dt      [20]Epot            [21]dEpot       [22]dEpot/dt
# [23]dErad         [24]dErad/dt        [25]E_bulk_BF   [26]dE_bulk_corr    [27] ([23]/dt)
# [28]E_bulk_AF     [29]dE_bulk_AF      [30] ([24]/dt)

mach = 10
clen = 1
ttc = clen / mach

fps = sys.argv[1::2] # files paths
lgs = sys.argv[2::2] # legends labels
    
# dissipation rate
for fp,lg in zip(fps,lgs):
    data = np.load(fp)

    t    = data[:,0] # time
    dK   = data[:,5] # dEkin(diss)

    i0 = (np.abs(t - 1*ttc)).argmin()

    xs = t[i0:]
    ys = -dK[i0:]

    xs = xs / ttc
    ys = np.cumsum(ys)

    lg += ' diss'
    plt.plot(xs, ys, '-', lw=3, label=lg)

# injection rate
plt.gca().set_prop_cycle(plt.matplotlib.rcParams['axes.prop_cycle'])
for fp,lg in zip(fps,lgs):
    data = np.load(fp)

    t  = data[:,0] # time
    dK = data[:,7] # dEkin(force)

    i0 = (np.abs(t - 1*ttc)).argmin()

    xs = t[i0:]
    ys = dK[i0:]

    xs = xs / ttc
    ys = np.cumsum(ys)

    lg += ' inject'
    plt.plot(xs, ys, '--', lw=3, label=lg)

plt.title("Turbulent Box (mach = %d): Evolution of Cummulated Kinetic Energy Dissipation" % mach)
plt.xlabel('characteristic time scale: t/t_c')
plt.ylabel('cummulative kinetic energy dissipation rate: Σ(t0,t)|ΔK|')

plt.xlim(0.8,3)
plt.ylim(0,230)
plt.grid()
plt.legend(loc='upper left')
#plt.legend(loc='upper right')

plt.show()

# fig = matplotlib.pyplot.gcf()
# fig.set_size_inches(18.5, 10.5)
# plt.tight_layout()
# plt.savefig(sys.stdout.buffer,format='png', bbox_inches='tight')
