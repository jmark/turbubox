#!/usr/bin/env python3

import sys
from matplotlib import pylab as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import imp
import gausslobatto
import flash, flexi, hopr
import flash_to_flexi as flfl
import scipy.misc
import ulz
import interpolate
import glob

sys.argv.reverse()
progpath = sys.argv.pop()
datafp = sys.argv.pop()
sinkfp = sys.argv.pop()

data = np.genfromtxt(datafp)

rs = data[:,0]
ps = data[:,1]

plt.grid()

xs = rs
ys = ps
plt.loglog(xs,ys, '-o', label='original')

xs = rs
ys = ps/rs**(-5/3)
plt.loglog(xs,ys, '-o', label='compensated: -5/3')

xs = rs
ys = ps/rs**(-2)
plt.loglog(xs,ys, '-o', label='compensated: -2')



plt.legend(loc='lower left')
plt.savefig(sinkfp,bbox_inches='tight')
