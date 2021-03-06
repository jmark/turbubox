import numpy as np
import ctypes as ct
from numpy.ctypeslib import ndpointer
import sys
import ulz

lib = ct.cdll.LoadLibrary(ulz.find_file('libshellavg.so', sys.path))

def shell_avg(X, nsamples=None, want_powerspectrum=False):
    if len(X.shape) == 2:
        return shell_avg_2d(X,nsamples,want_powerspectrum)
    if len(X.shape) == 3:
        return shell_avg_3d(X,nsamples,want_powerspectrum)
    raise NotImplementedError("{} dimensions is not supported.".format(len(X.shape)))

lib.shell_avg_2d.argtypes = [
    # void shell_avg_3d(
    #     const double *X, const int Nx, const int Ny,
    #     double *rs, double *cs, double *ts, const int nsamples
    # );

    ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ct.c_int, ct.c_int,

    ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ct.c_int, ct.c_int
]
 
def shell_avg_2d(X, nsamples=None,want_powerspectrum=False):
    Nx,Ny = X.shape

    if not nsamples: nsamples = min(Nx,Ny)

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
    lib.shell_avg_2d(X,Nx,Ny, cs,rs,ts,nsamples,int(want_powerspectrum))

    # take average and return
    # note: if div-by-zero warning arises: the inputs-to-samples ratio 
    # is not adequate

    if want_powerspectrum:
        return rs/cs, 2*np.pi*ts/cs
    return rs/cs, ts/cs

lib.shell_avg_2d_corner.argtypes = [
    # void shell_avg_3d(
    #     const double *X, const int Nx, const int Ny,
    #     double *rs, double *cs, double *ts, const int nsamples
    # );

    ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ct.c_int, ct.c_int,

    ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ct.c_int, ct.c_int
]
 
def shell_avg_2d_corner(X, nsamples=None,want_powerspectrum=False):
    Nx,Ny = X.shape

    if not nsamples: nsamples = min(Nx,Ny)

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
    lib.shell_avg_2d_corner(X,Nx,Ny, cs,rs,ts,nsamples,int(want_powerspectrum))

    # take average and return
    # note: if div-by-zero warning arises: the inputs-to-samples ratio 
    # is not adequate

    if want_powerspectrum:
        return rs/cs, 2*np.pi*ts/cs
    return rs/cs, ts/cs

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
    ct.c_int, ct.c_int
]
 
def shell_avg_3d(X, nsamples=None,want_powerspectrum=False):
    Nx,Ny,Nz = X.shape

    if not nsamples: nsamples = min(Nx,Ny,Nz)

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
    lib.shell_avg_3d(X,Nx,Ny,Nz, cs,rs,ts,nsamples,int(want_powerspectrum))

    # take average and return
    # note: if div-by-zero warning arises: the inputs-to-samples ratio 
    # is not adequate

    if want_powerspectrum:
        return rs/cs, 4*np.pi*ts/cs
    return rs/cs, ts/cs

lib.shell_avg_3d_corner.argtypes = [
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ct.c_int, ct.c_int, ct.c_int,

    ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ct.c_int, ct.c_int
]
 
def shell_avg_3d_corner(X, nsamples=None,want_powerspectrum=False):
    Nx,Ny,Nz = X.shape

    if not nsamples: nsamples = min(Nx,Ny,Nz)

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
    lib.shell_avg_3d_corner(X,Nx,Ny,Nz, cs,rs,ts,nsamples,int(want_powerspectrum))

    # take average and return
    # note: if div-by-zero warning arises: the inputs-to-samples ratio 
    # is not adequate

    if want_powerspectrum:
        return rs/cs, 4*np.pi*ts/cs
    return rs/cs, ts/cs



lib.shell_max_3d.argtypes = [
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
 
def shell_max_3d(X, nsamples=None):
    Nx,Ny,Nz = X.shape

    if not nsamples: nsamples = min(Nx,Ny,Nz)

    # setup result arrays
    cs = np.zeros(nsamples, dtype=np.double) # counts
    rs = np.zeros(nsamples, dtype=np.double) # radii
    #ts = np.full(nsamples, np.finfo(np.double).min) # totals
    ts = np.full(nsamples, -99) # totals

    # ensure C-like arrays
    X  = np.require( X.ravel(), dtype=np.double, requirements=['C', 'A'])
    cs = np.require(cs.ravel(), dtype=np.double, requirements=['C', 'A', 'W'])
    rs = np.require(rs.ravel(), dtype=np.double, requirements=['C', 'A', 'W'])
    ts = np.require(ts.ravel(), dtype=np.double, requirements=['C', 'A', 'W'])

    # call C-function
    lib.shell_max_3d(X,Nx,Ny,Nz, cs,rs,ts,nsamples)

    # take average and return
    # note: if div-by-zero warning arises: the inputs-to-samples ratio 
    # is not adequate

    return rs/cs, ts

if __name__ == '__main__':
    #X = np.random.rand(100,100,100)
    #X = np.ones((100,100,100))
    X = np.ones((100,100))
    
    #rs,ts = shell_avg_3d(X)
    rs,ts = shell_avg_2d(X)

    import sys
    np.savetxt(sys.stdout.buffer, np.array([rs,ts,rs**2 * ts]).T)
