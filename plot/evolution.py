#!/usr/bin/env python3

import sys
import pickle
import numpy as np
from matplotlib import pylab as plt

fps = sys.argv[1::2]
lgs = sys.argv[2::2]

keys = 'taskid step time turn mach ekintot vorttot'.split()
b3fp = '../../data/evol-es.pickle'

with open(b3fp, 'rb') as fd:
    data = pickle.load(fd)

data = data[:194].T
xs = data[2]
ys = data[5]

n = 2
b3ys = np.abs((np.roll(ys,n) - np.roll(ys,-n)) / (np.roll(xs,n) - np.roll(xs,-n)))

for fp,lg in zip(fps,lgs):
    with open(fp, 'rb') as fd:
        data = pickle.load(fd)

    data = data[:194].T
    xs = data[2]
    ys = data[5]

    n = 2
    ys = np.abs((np.roll(ys,n) - np.roll(ys,-n)) / (np.roll(xs,n) - np.roll(xs,-n)))

    ys = (ys / b3ys)

    #xs = xs[80:-1]
    xs = data[4][80:-1]
    ys = ys[80:-1]

    #xs = list(reversed(xs))
    #ys = list(reversed(ys))

    #print(list(reversed(xs)))
    #break 
    #xs = np.log10(xs)
    #ys = np.log10(ys)

    plt.plot(xs,ys, '-', lw=2, label=lg)

plt.grid()
plt.legend()
plt.show()
