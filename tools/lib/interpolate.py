import numpy as np
import ctypes as ct
from numpy.ctypeslib import ndpointer
import sys
import ulz

lib = ct.cdll.LoadLibrary(ulz.find_file('libinterpolate.so', sys.path))

# =========================================================================== #
# double
# LagrangePolynomial(const double *xs, const int xslen, const int j, const double X);

lib.LagrangePolynomial.restype = ct.c_double
lib.LagrangePolynomial.argtypes = [
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"), ct.c_int, ct.c_int, ct.c_double
]

def LagrangePolynomial(xs, j, X):
    xs = np.require(xs.ravel(), dtype=np.double, requirements=['C', 'A'])
    return lib.LagrangePolynomial(xs, len(xs), int(j), float(X))

# =========================================================================== #

# void
# lagrange_interpolate_2d_RG(
#     const int xslen, const double *xs, const double *fs,
#     const int Xslen, const double *Xs,       double *Fs
# );

lib.lagrange_interpolate_2d_RG.argtypes = [
    ct.c_int, ndpointer(ct.c_double, flags="C_CONTIGUOUS"), ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ct.c_int, ndpointer(ct.c_double, flags="C_CONTIGUOUS"), ndpointer(ct.c_double, flags="C_CONTIGUOUS")
]

def lagrange_interpolate_2d_RG(xs, Xs, fs):
    shape = fs.shape

    xs = np.require(xs.ravel(), dtype=np.double, requirements=['C','A'])
    Xs = np.require(Xs.ravel(), dtype=np.double, requirements=['C','A'])
    fs = np.require(fs.ravel(), dtype=np.double, requirements=['C','A'])

    Fs = np.zeros(len(Xs)**2,dtype=np.double)
    Fs = np.require(Fs.ravel(), dtype=np.double, requirements=['C','A','W'])

    lib.lagrange_interpolate_2d_RG(len(xs),xs,fs, len(Xs),Xs,Fs)

    return Fs.reshape([len(Xs)]*2)

# =========================================================================== #

# void
# lagrange_interpolate_3d_RG(
#     const int xslen, const double *xs, const double *fs,
#     const int Xslen, const double *Xs,       double *Fs
# );

lib.lagrange_interpolate_3d_RG.argtypes = [
    ct.c_int, ndpointer(ct.c_double, flags="C_CONTIGUOUS"), ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ct.c_int, ndpointer(ct.c_double, flags="C_CONTIGUOUS"), ndpointer(ct.c_double, flags="C_CONTIGUOUS")
]

def lagrange_interpolate_3d_RG(xs, Xs, fs):
    shape = fs.shape

    xs = np.require(xs.ravel(), dtype=np.double, requirements=['C','A'])
    Xs = np.require(Xs.ravel(), dtype=np.double, requirements=['C','A'])
    fs = np.require(fs.ravel(), dtype=np.double, requirements=['C','A'])

    Fs = np.zeros(len(Xs)**3,dtype=np.double)
    Fs = np.require(Fs.ravel(), dtype=np.double, requirements=['C','A','W'])

    lib.lagrange_interpolate_3d_RG(len(xs),xs,fs, len(Xs),Xs,Fs)

    return Fs.reshape([len(Xs)]*3)

# =========================================================================== #

# void
# flash_to_flexi_RG(
#     const int xslen, const double *xs, const double *fss,
#     const int Xslen, const double *Xs,       double *Fss,
#     const int oflen, const int *offsets
# );

# flash_to_flexi_RG(
#     const int xslen, const double *xs, 
#     const int Nx, const int Ny, const int Nz, const double *fss,
#     const int Xslen, const double *Xs,       double *Fss,
#     const int oflen, const int *offsets
# )

lib.box_to_flexi.argtypes = [
    ct.c_int, ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ct.c_int, ct.c_int, ct.c_int, ndpointer(ct.c_double, flags="C_CONTIGUOUS"), 
    ct.c_int, ndpointer(ct.c_double, flags="C_CONTIGUOUS"), ndpointer(ct.c_double, flags="C_CONTIGUOUS"), 
    ct.c_int, ndpointer(ct.c_int, flags="C_CONTIGUOUS")
]

def box_to_flexi(xs, Xs, box, flx):
    shape = box.shape

    # ll ... lower left
    # tr ... top right
    lls, trs = flx.mesh.get_cell_coords()
    Is,Js,Ks = tuple(np.round((flx.npoly+1) * lls / flx.mesh.cellsize).astype(int).T)
    offsets  = ((Is * shape[1]) + Js) * shape[2] + Ks
    offsets  = np.require(offsets.ravel(), dtype=np.int32, requirements=['C', 'A'])

    xs  = np.require(xs.ravel(), dtype=np.double, requirements=['C', 'A'])
    Xs  = np.require(Xs.ravel(), dtype=np.double, requirements=['C', 'A'])
    box = np.require(box.ravel(), dtype=np.double, requirements=['C', 'A'])

    flxdata = np.empty(len(offsets) * len(Xs)**3, dtype=np.double)
    flxdata = np.require(flxdata.ravel(), dtype=np.double, requirements=['C', 'A', 'W'])

    lib.box_to_flexi(
        len(xs), xs,
        *shape, box,
        len(Xs), Xs, flxdata,
        len(offsets), offsets)

    return flxdata.reshape(len(offsets),*[len(Xs)]*3)

# =========================================================================== #

# void
# flexi_to_box(
#     const int xslen, const double *xs,
#     const int Xslen, const double *Xs, 
#     const int nelems, const int *offsets, double *flexi
#     const int Nx, const int Ny, const int Nz, const double *box,
# );

lib.flexi_to_box.argtypes = [
    ct.c_int, ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ct.c_int, ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ct.c_int, ndpointer(ct.c_int, flags="C_CONTIGUOUS"), ndpointer(ct.c_double, flags="C_CONTIGUOUS"), 
    ct.c_int, ct.c_int, ct.c_int, ndpointer(ct.c_double, flags="C_CONTIGUOUS")
]

def flexi_to_box(xs, Xs, flxdata, flx):
    shape = tuple(len(Xs) * flx.mesh.gridsize.astype(np.int))

    # ll ... lower left
    # tr ... top right
    lls, trs = flx.mesh.get_cell_coords()
    Is,Js,Ks = tuple(np.round(sh*ll) for sh,ll in zip(shape,lls.T))
    offsets  = ((Is * shape[1]) + Js) * shape[2] + Ks
    offsets  = np.require(offsets.ravel(), dtype=np.int32, requirements=['C', 'A'])

    xs = np.require(xs.ravel(), dtype=np.double, requirements=['C', 'A'])
    Xs = np.require(Xs.ravel(), dtype=np.double, requirements=['C', 'A'])

    flxdata = flxdata.transpose(0,3,2,1)
    flxdata = np.require(flxdata.ravel(), dtype=np.double, requirements=['C', 'A'])

    box = np.zeros(shape,dtype=np.double)
    box = np.require(box.ravel(), dtype=np.double, requirements=['C', 'A', 'W'])

    lib.flexi_to_box(
        len(xs), xs,
        len(Xs), Xs,
        len(offsets), offsets, flxdata,
        *shape, box)

    return box.reshape(shape)
