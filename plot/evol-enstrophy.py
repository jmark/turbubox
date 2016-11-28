#!/usr/bin/env python3

import sys
import pickle
import numpy as np
from matplotlib import pylab as plt
import matplotlib
matplotlib.rc('font', family='DejaVu Sans')
matplotlib.rcParams.update({'font.size': 20})

import cycler
ccycle = cycler.cycler('color', ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00'])

plt.rc('axes', prop_cycle=ccycle)

# forcing file
fps = sys.argv[0:7][1::2]
lgs = sys.argv[0:7][2::2]

mach = 10
clen = 1
ttc = clen / mach

for fp,lg in zip(fps,lgs):
    with open(fp, 'rb') as fd:
        data = pickle.load(fd)

    data = data.T

    t = data[2]
    E = data[6]

    xs = t / ttc
    ys = E
    #lg += ' (snapshots)'
    plt.plot(xs,ys, '-', lw=3, label=lg)

plt.title("Turbulent Box (Mach = %d): Evolution of Total Enstrophy" % mach)
plt.xlabel('characteristic time scale: t/t_c')
plt.ylabel('Total Enstrophy E')

plt.xlim(0,5)

plt.grid()
plt.legend()

plt.show()

# fig = matplotlib.pyplot.gcf()
# fig.set_size_inches(18.5, 10.5)
# plt.tight_layout()
# plt.savefig(sys.stdout.buffer,format='png', bbox_inches='tight')
