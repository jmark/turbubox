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

def mymean(xs):
    avg = np.mean(xs)
    return np.mean([x for x in xs if np.abs((avg-x)/avg) <= 0.5])

def moving_avg(xs, n=2):
    return np.array([mymean(x) for x in xs[:xs.shape[0]//n * n].reshape(-1,n)])


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

dname = 'vort'
start = 30
end = 80

from scipy.optimize import curve_fit

def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x-mu)**2/2/sigma)

def lognormal(x, A, mu, sigma):
    return A * np.exp(mu + sigma * (x-B))

def loglognormal(x, A, mu, sigma):
    return np.log10(A * np.exp(mu + sigma * np.log10(x)))

for fp,lg in zip(fps,lgs):
    with open(fp, 'rb') as fd:
        data = pickle.load(fd)

    rs, ps = data[start][dname]
    ps /= np.sum(ps)

    rsg, psg = np.zeros_like(rs), np.zeros_like(ps)

    cnt = 1
    for x in data[start+11:end]:
        cnt += 1
        r,p = x[dname]
        rs += r
        ps += p

    rs /= cnt
    ps /= cnt

    cnt = 1
    for x in data[start+11:end]:
        cnt += 1
        r,p = x[dname]
        rsg += (r - rs)**2
        psg += (p - ps)**2

    xs = ps[:-1] + (ps[1]-ps[0])/2
    ys = rs

    sigma = np.sqrt(rsg/(cnt-1))
    ysT = (ys + sigma)
    ysD = (ys - sigma)

    n = 20
    xsB = moving_avg( xs, 5)[:-n]
    ysT = moving_avg(ysT, 5)[:-n]
    ysD = moving_avg(ysD, 5)[:-n]

    plt.fill_between(xsB, ysD, ysT, color='gray', alpha=0.3)
    #plt.semilogy(xs,ys, '-', lw=3, label=lg)
    plt.plot(xs,ys, '-', lw=3, label=lg)

# plt.gca().set_prop_cycle(plt.matplotlib.rcParams['axes.prop_cycle'])
# 
# for fp,lg in zip(fps,lgs):
#     with open(fp, 'rb') as fd:
#         data = pickle.load(fd)
# 
#     rs, ps = data[start][dname]
#     ps /= np.sum(ps)
# 
#     rsg, psg = np.zeros_like(rs), np.zeros_like(ps)
# 
#     cnt = 1
#     for x in data[start+11:end]:
#         cnt += 1
#         r,p = x[dname]
#         rs += r
#         ps += p
# 
#     rs /= cnt
#     ps /= cnt
# 
#     xs = ps[:-1] + (ps[1]-ps[0])/2
#     ys = rs
# 
#     popt, pcov = curve_fit(gaussian, xs, ys)
#     ysF = gaussian(xs, *popt)
# 
#     plt.semilogy(xs,ysF, '--', lw=3, label=lg + ' gaussian fit: A=%4.2f, mu=%4.2f, sigma^2=%4.2f' % tuple(popt))

plt.title("Turbulent Box (Mach = %d): Velocity PDF" % mach)
plt.xlabel('velocity')
plt.ylabel('velocity pdf')

plt.grid()
plt.legend(loc='lower center')

plt.show()

# fig = matplotlib.pyplot.gcf()
# fig.set_size_inches(18.5, 10.5)
# plt.tight_layout()
# plt.savefig(sys.stdout.buffer,format='png', bbox_inches='tight')
