#!/usr/bin/env python3

import sys
import pickle
import numpy as np
from matplotlib import pylab as plt
import matplotlib
matplotlib.rc('font', family='DejaVu Sans')
matplotlib.rcParams.update({'font.size': 20})

# [01]time           [02]dt  [03]EkinLastStp   [04]EkinBefFor   [05]EkinAftFor
# [06]dEkin(diss)   [07] ([06]/dt)  [08]dEkin(forc)   [09] ([08]/dt)
# [10]Vrms_vw_BF   [11]Vrms_vw_AF   [12]Vrms_mw_BF   [13]Vrms_mw_AF
# [14]Eint        [15]dEint     [16]dEint/dt         [17]Emag        [18]dEmag
# [19]dEmag/dt         [20]Epot        [21]dEpot     [22]dEpot/dt
# [23]dErad     [24]dErad/dt    [25]E_bulk_BF [26]dE_bulk_corr   [27] ([23]/dt)
# [28]E_bulk_AF   [29]dE_bulk_AF   [30] ([24]/dt)

start = {'b3': 464, 'b5': 501, 'es': 9335}

mach = 10
clen = 1
ttc = clen / mach

fps = sys.argv[1::2] # files paths
lgs = sys.argv[2::2] # legends labels

for fp,lg in zip(fps[:-1],lgs[:-1]):
    data = np.load(fp)
    data = data[start[lg]:]

    t  = data[:,0] # time
    dt = data[:,1] # time step
    dkd = data[:,5] # dEkin(diss)
    dkf = data[:,7] # dEkin(force)

    xs = t / ttc

    ys = dkd
    ys = np.cumsum(ys) 
    ys = np.abs(ys)
    plt.plot(xs, ys, '-', lw=3, label=lg + ' dissip')

    ys = dkf
    ys = np.cumsum(ys) 
    plt.plot(xs, ys, '-', lw=3, label=lg + ' inject')


fp = sys.argv[-2]
lg = sys.argv[-1]

data = np.load(fp)
data = data[start[lg]:]

t  = data[:,0] # time
dt = data[:,1] # time step
dkd = data[:,5] # dEkin(diss)
dkf = data[:,7] # dEkin(force)

xs = t / ttc

ys = dkd
ys = np.cumsum(ys) 
ys = np.abs(ys)
plt.plot(xs, ys, '-', lw=3, label=lg + ' dissip')

ys = dkf
ys = np.cumsum(ys) 
plt.plot(xs, ys, '-', lw=3, label=lg + ' inject')

plt.title("Turbulent Box (mach = %d): Evolution of Cummulative Dissipated and Injected Kinetic Energy" % mach)
plt.xlabel('characteristic time scale: t/t_c')
plt.ylabel('Cummulative Kinetic Energy Dissipation/Injection: Σ|ΔK/Δt|(t)')

plt.xlim(0.75,3)
plt.ylim(-5,250)
plt.grid()
plt.legend(loc='upper left')
#plt.legend(loc='upper right')

plt.show()

# fig = matplotlib.pyplot.gcf()
# fig.set_size_inches(18.5, 10.5)
# plt.tight_layout()
# plt.savefig(sys.stdout.buffer,format='png', bbox_inches='tight')
