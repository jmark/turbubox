#!/usr/bin/env python3

import sys
import pickle
import numpy as np
from matplotlib import pylab as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 20})

# [01]time           [02]dt  [03]EkinLastStp   [04]EkinBefFor   [05]EkinAftFor
# [06]dEkin(diss)   [07] ([06]/dt)  [08]dEkin(forc)   [09] ([08]/dt)
# [10]Vrms_vw_BF   [11]Vrms_vw_AF   [12]Vrms_mw_BF   [13]Vrms_mw_AF
# [14]Eint        [15]dEint     [16]dEint/dt         [17]Emag        [18]dEmag
# [19]dEmag/dt         [20]Epot        [21]dEpot     [22]dEpot/dt
# [23]dErad     [24]dErad/dt    [25]E_bulk_BF [26]dE_bulk_corr   [27] ([23]/dt)
# [28]E_bulk_AF   [29]dE_bulk_AF   [30] ([24]/dt)

mach = 10
clen = 1
ttc = clen / mach

fps = sys.argv[1::2] # files paths
lgs = sys.argv[2::2] # legends labels

for fp,lg in zip(fps,lgs):
    data = np.load(fp)

    time = data[:,0] # time
    dt   = data[:,1] # time
    dk   = data[:,7] # dEkinForce

    #xs = (time + dt/2)[:-1]
    #ys = (dt * (dk + np.roll(dk,-1))/2)[:-1]
    #ys = np.cumsum(ys) 

    xs = time / ttc
    ys = dt / ttc

    plt.plot(xs, ys, '-', lw=3, label=lg)

fp = sys.argv[-2]
lg = 'es: lifted by 10*dt_es'
data = np.load(fp)

time = data[:,0] # time
dt   = data[:,1] # time
dk   = data[:,7] # dEkinForce

#xs = (time + dt/2)[:-1]
#ys = (dt * (dk + np.roll(dk,-1))/2)[:-1]
#ys = np.cumsum(ys) 

xs = time / ttc
ys = 10 * dt / ttc

plt.plot(xs, ys, '-', lw=4, label=lg)

plt.title("Turbulent Box (mach = %d): Evolution of Time Step 'dt'" % mach)

plt.xlabel('characteristic time scale: t/t_c')
plt.ylabel('characteristic time step scale: dt/t_c')

plt.xlim(0,5)
plt.ylim(0,0.011)
plt.grid()
#plt.legend(loc='upper left')
plt.legend(loc='upper right')
#plt.show()

fig = matplotlib.pyplot.gcf()
fig.set_size_inches(18.5, 10.5)
plt.tight_layout()
plt.savefig(sys.stdout.buffer,format='png', bbox_inches='tight')
