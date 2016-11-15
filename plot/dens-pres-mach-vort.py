#!/usr/bin/env python3

import sys
import numpy as np
import ulz
import flash

sys.argv.reverse()
progname = sys.argv.pop()
flashfp  = sys.argv.pop()
ifile    = int(sys.argv.pop())
nfile    = int(sys.argv.pop())
sinkfp   = sys.argv.pop()

fls = flash.File(flashfp)

time = fls.realscalars['time']
step = fls.integerscalars['nstep']
c_s  = fls.realruntime['c_ambient']
rho0 = fls.realruntime['rho_ambient']
LEN  = fls.domainsize[0]

turntime = time / (LEN / c_s / 10)

GS = fls.gridsize
DS = fls.domainsize
CS = fls.cellsize

Vgrid   = np.prod(GS)
Vcell   = np.prod(CS) 
Vdomain = np.prod(DS) 

dens = fls.data('dens')
pres = fls.data('pres')
vels = [fls.data('vel'+dim) for dim in 'x y z'.split()]
mach = np.sqrt(ulz.norm(*vels)/3) / c_s
vort = CS[0]**5/12.0 * dens * ulz.norm(*ulz.curl(vels[0],vels[1],vels[2],CS[0],CS[1],CS[2]))

ax = 2
cdens = np.log10(np.sum(dens,axis=ax))
cpres = np.log10(np.sum(pres,axis=ax))
cmach = np.sum(mach,axis=ax)/mach.shape[ax]
cvort = np.log10(np.sum(vort,axis=ax))

# print(ifile, 
#     np.max(cdens),
#     np.max(cpres),
#     np.max(cmach),
#     np.max(cvort))

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

matplotlib.rcParams.update({'font.size': 10})

subplt = [2,2,0]
fig = plt.figure(figsize=(15, 12))

st = plt.suptitle(
    "Stirred turbulence in periodic box: mach 10 | B5 | t_c = %1.3f (frame: %03d/%03d)" % (turntime, ifile, nfile),
    fontsize='x-large')
st.set_y(1.01)

def plot(data, title, crange=None):
    subplt[2] += 1
    ax = fig.add_subplot(*subplt)
    ax.set_title('column %s' % title)
    ax.set_xlabel('x index')
    ax.set_ylabel('y index')
    if crange: img = ax.imshow(data, cmap=plt.get_cmap('cubehelix'),**crange)
    else: img = ax.imshow(data, cmap=plt.get_cmap('cubehelix'))
    plt.colorbar(img,fraction=0.046, pad=0.04, format='%1.2f')

crange = {'vmin': 0, 'vmax': 4}
plot(cdens, 'density (log10)', crange)

crange = {'vmin': 0, 'vmax': 4}
plot(cpres, 'pressure (log10)', crange)

crange = {'vmin': 0, 'vmax': 10}
plot(cmach, 'sonic mach number (grid normalized)', crange)

crange = {'vmin': -10, 'vmax': -3.4}
plot(cvort, 'vorticity (log10)', crange)

fig.tight_layout()
plt.savefig(sinkfp,bbox_inches='tight')
