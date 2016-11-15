#!/usr/bin/env python3

import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import ulz
import flash

matplotlib.rcParams.update({'font.size': 10})

sys.argv.reverse()
progname = sys.argv.pop()
flashfpB3 = sys.argv.pop()
flashfpB5 = sys.argv.pop()
ifile  = int(sys.argv.pop())
nfile  = int(sys.argv.pop())
sinkfp = sys.argv.pop()

fls = flash.File(flashfpB3)

time = fls.realscalars['time']
step = fls.integerscalars['nstep']

c_s  = fls.realruntime['c_ambient']
rho0 = fls.realruntime['rho_ambient']

LEN  = fls.domainsize[0]

turntime = time / (LEN / c_s / 10)

subplt = [1,2,0]
crange = {'vmin': 0, 'vmax': 4}
#crange = {'vmin': 0, 'vmax': 4000}

fig = plt.figure(figsize=(15, 12))

st = plt.suptitle("Stirred turbulence in periodic box: mach 10", fontsize='x-large')
st.set_y(0.8)

def plot(data, title, crange=None):
    subplt[2] += 1
    ys = np.log10(np.sum(data,axis=2)).T
    #ys = np.sum(data ,axis=2).T
    ax = fig.add_subplot(*subplt)
    ax.set_title('column density: %s / t_c = %1.2f (frame: %03d/%03d)' % title)
    ax.set_xlabel('x index')
    ax.set_ylabel('y index')
    if crange: img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'),**crange)
    else: img = ax.imshow(ys, cmap=plt.get_cmap('cubehelix'))
    plt.colorbar(img,fraction=0.046, pad=0.04, format='%1.2f')

fls = flash.File(flashfpB3)
plot(fls.data('dens'), ('B3', turntime, ifile, nfile), crange)

fls = flash.File(flashfpB5)
plot(fls.data('dens'), ('B5', turntime, ifile, nfile), crange)

fig.tight_layout()
plt.savefig(sinkfp,bbox_inches='tight')
