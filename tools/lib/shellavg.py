import numpy as np
import ctypes as ct
from numpy.ctypeslib import ndpointer
import sys
import ulz

lib = ct.cdll.LoadLibrary(ulz.find_file('libshellavg.so', sys.path))

lib.shell_avg_3d.argtypes = [
    # void shell_avg_3d(
    #     const double *X, const int Nx, const int Ny, const int Nz,
    #     double *rs, double *cs, double *ts, const int nsamples
    # );

    ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ct.c_int, ct.c_int, ct.c_int,

    ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ct.c_int
]
 
def shell_avg_3d(X, nsamples=None):
    Nx,Ny,Nz = X.shape

    if not nsamples:
        nsamples = min(Nx,Ny,Nz)

    # setup result arrays
    cs = np.zeros(nsamples, dtype=np.double) # counts
    rs = np.zeros(nsamples, dtype=np.double) # radii
    ts = np.zeros(nsamples, dtype=np.double) # totals

    # ensure C-like arrays
    X  = np.require( X.ravel(), dtype=np.double, requirements=['C', 'A'])
    cs = np.require(cs.ravel(), dtype=np.double, requirements=['C', 'A', 'W'])
    rs = np.require(rs.ravel(), dtype=np.double, requirements=['C', 'A', 'W'])
    ts = np.require(ts.ravel(), dtype=np.double, requirements=['C', 'A', 'W'])

    # call C-function
    lib.shell_avg_3d(X,Nx,Ny,Nz, cs,rs,ts,nsamples)

    # take average
    # note: if div-by-zero warning arises: the inputs-to-samples ratio 
    # is not adequate
    rs /= cs
    ts /= cs

    return rs,ts

if __name__ == '__main__':
    #X = np.random.rand(100,100,100)
    X = np.ones((100,100,100))
    
    rs,ts = shell_avg_3d(X)

    import sys
    np.savetxt(sys.stdout.buffer, np.array([rs,ts,rs**2 * ts]).T)
