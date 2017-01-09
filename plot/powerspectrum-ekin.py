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

fps = sys.argv[0:7][1::2]
lgs = sys.argv[0:7][2::2]

mach = 10
clen = 1
ttc = clen / mach

dname = 'ekin'
start = 40
end = 80

for fp,lg in zip(fps,lgs):
    with open(fp, 'rb') as fd:
        data = pickle.load(fd)

    rs, ps = data[start][dname]
    rsg, psg = np.zeros_like(rs), np.zeros_like(ps)

    cnt = 1
    for x in data[start:end]:
        cnt += 1
        r,p = x[dname]
        rs += r
        ps += p

    rs /= cnt
    ps /= cnt

    cnt = 1
    for x in data[start:end]:
        cnt += 1
        r,p = x[dname]
        rsg += (r - rs)**2
        psg += (p - ps)**2

    rsg = np.sqrt(rsg/(cnt-1))
    psg = np.sqrt(psg/(cnt-1))

    xs = rs[1:]
    ys = xs**2 * ps[1:]

    ysT = xs**2 * (ps + psg)[1:]
    ysD = xs**2 * (ps - psg)[1:]

    plt.fill_between(xs, ysD, ysT, color='gray', alpha=0.3)

    plt.loglog(xs,ys, '-', lw=3, label=lg)

plt.gca().set_prop_cycle(plt.matplotlib.rcParams['axes.prop_cycle'])
for fp,lg in zip(fps,lgs):
    with open(fp, 'rb') as fd:
        data = pickle.load(fd)

    rs, ps = data[start][dname]
    cnt = 1
    for x in data[start+11:end]:
        cnt += 1
        r,p = x[dname]
        rs += r
        ps += p

    rs /= cnt
    ps /= cnt

    corr = -0.62
    xs = rs[1:]
    ys = xs**2 * ps[1:] / xs**corr

    lg += ': P(K)/k^%3.2f' % corr
    plt.loglog(xs,ys, '--', lw=3, label=lg)

plt.title("Turbulent Box (Mach = %d): Total Kinetic Energy Power Spectrum" % mach)
plt.xlabel('grid number: k')
plt.ylabel('kinetic energy spectrum: P(K)')

plt.grid()
plt.legend(loc='lower left')

plt.show()

# fig = matplotlib.pyplot.gcf()
# fig.set_size_inches(18.5, 10.5)
# plt.tight_layout()
# plt.savefig(sys.stdout.buffer,format='png', bbox_inches='tight')
