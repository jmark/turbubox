#!/usr/bin/env python3

import sys
from matplotlib import pylab as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import gausslobatto as gl
import flash, flexi, hopr
import flash_to_flexi as flfl
import scipy.misc
import ulz

sys.argv.reverse()
sys.argv.pop()
outpath = sys.argv.pop()

lin1d = np.linspace(0.1,0.9,4)
func  = lambda x: np.sin(2*np.pi*x)

Lv = lambda x: np.array([gl.LagrangePolynomial(lin1d,j,x) for j in range(len(lin1d))])
Lf = lambda f,x: np.dot(f,Lv(x))

plt.figure(figsize=(10, 5))
xs = np.linspace(-0.2,1.2,30)
plt.grid()
plt.plot(xs, func(xs), '-o', lw=2)
plt.plot(xs, [Lf(func(lin1d),x) for x in xs], '-o', lw=2)
plt.plot(lin1d, func(lin1d), 'o', lw=2, markersize=10)
plt.legend(['original: sin(2 Pi x)', 'interpolated', 'sample nodes'], loc='upper left')

plt.savefig(outpath,bbox_inches='tight')
