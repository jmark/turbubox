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

t_ndouble = ndpointer(ct.c_double, flags="C_CONTIGUOUS")
t_int = ct.c_int

# void
# box_to_elements(
#     int Nx, int Ny, int Nz, double *boxptr, 
#     int nelems, int nx, int ny, int nz, double *offsetsptr, double *elemsptr);

lib.box_to_elements.argtypes = [
    t_int, t_int, t_int, t_ndouble,
    t_int, t_ndouble, 
    t_int, t_int, t_int, t_ndouble, t_int
]

def box_to_elements(box, flx, neighbors=0):
    Nnodes   = (np.array(box.shape)//flx.mesh.meshshape)[0]

    lls, _   = flx.mesh.get_cell_coords()
    offsets  = Nnodes * lls / flx.mesh.elemsize

    N = Nnodes + 2*neighbors
    elems    = np.zeros([flx.mesh.nrelems, N,N,N], dtype=np.double)

    boxptr   = np.require(box.ravel(), dtype=np.double, requirements=['C','A'])
    elemsptr = np.require(elems.ravel(), dtype=np.double, requirements=['C','A'])
    ofsptr   = np.require(offsets.ravel(), dtype=np.double, requirements=['C', 'A'])

    lib.box_to_elements(*box.shape, boxptr, flx.mesh.nrelems, ofsptr, *elems[0].shape, elemsptr, neighbors)

    return elemsptr.reshape(elems.shape)

# =========================================================================== #

lib.elements_to_box.argtypes = [
    t_int, t_int, t_int, t_ndouble,
    t_int, t_ndouble, 
    t_int, t_int, t_int, t_ndouble
]

def elements_to_box(elems, mesh):
    # lower left corners normed to unit intervall
    lls      = (mesh.domain[0] + mesh.elemcoords[0])/mesh.domainsize

    box      = np.zeros(elems[0].shape * mesh.meshshape)
    offsets  = np.array(box.shape) * lls

    boxptr   = np.require(box.ravel(), dtype=np.double, requirements=['C','A'])
    elemsptr = np.require(elems.ravel(), dtype=np.double, requirements=['C','A'])
    ofsptr   = np.require(offsets.ravel(), dtype=np.double, requirements=['C', 'A'])

    lib.elements_to_box(*box.shape, boxptr, mesh.nrelems, ofsptr, *elems[0].shape, elemsptr)

    return boxptr.reshape(box.shape)

# =========================================================================== #

lib.box_to_elements_avg_boundaries.argtypes = [
    t_int, t_int, t_int, t_ndouble,
    t_int, t_ndouble, 
    t_int, t_int, t_int, t_ndouble
]

def box_to_elements_avg_boundaries(box, flx):
    lls, _   = flx.mesh.get_cell_coords()
    offsets  = (flx.Nout) * lls / flx.mesh.cellsize

    N = flx.Nout
    elems    = np.zeros([flx.mesh.nrelems, N,N,N], dtype=np.double)

    boxptr   = np.require(box.ravel(), dtype=np.double, requirements=['C','A'])
    ofsptr   = np.require(offsets.ravel(), dtype=np.double, requirements=['C', 'A'])
    elemsptr = np.require(elems.ravel(), dtype=np.double, requirements=['C','A','W'])

    lib.box_to_elements_avg_boundaries(*box.shape, boxptr, flx.mesh.nrelems, ofsptr, *elems[0].shape, elemsptr)

    return elemsptr.reshape(elems.shape)

# =========================================================================== #

# void
# change_grid_space(
#     const int nelems,
#     const int nx, const int ny, const int nz ,const double *xs, double *fss, 
#     const int Nx, const int Ny, const int Nz ,const double *Xs, double *Fss);

lib.change_grid_space.argtypes = [
   t_int,
   t_int, t_int, t_int, t_ndouble, t_ndouble,
   t_int, t_int, t_int, t_ndouble, t_ndouble
]

def change_grid_space(fs,xs,Xs):
    Fs = np.empty([len(fs), len(Xs), len(Xs), len(Xs)])

    xsptr = np.require(xs, dtype=np.double, requirements=['C','A'])
    Xsptr = np.require(Xs, dtype=np.double, requirements=['C','A'])
    fsptr = np.require(fs, dtype=np.double, requirements=['C','A'])
    Fsptr = np.require(Fs, dtype=np.double, requirements=['C','A','W'])

    lib.change_grid_space(
        len(fs),
        len(xs), len(xs), len(xs), xsptr, fsptr,
        len(Xs), len(Xs), len(Xs), Xsptr, Fsptr
    )

    return Fsptr.reshape(Fs.shape)

# =========================================================================== #

# deprecated ...
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
    Is,Js,Ks = tuple(np.round((flx.Nout) * lls / flx.mesh.cellsize).astype(int).T)
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
# box_to_flexi_with_averaged_boundaries(
#     const int xslen, const double *xs, 
#     const int Nx, const int Ny, const int Nz, const double *fss,
#     const int Xslen, const double *Xs, double *Fss,
#     const int oflen, const int *offsets
# )

# lib.box_to_flexi_with_averaged_boundaries.argtypes = [
#     ct.c_int, ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
#     ct.c_int, ct.c_int, ct.c_int, ndpointer(ct.c_double, flags="C_CONTIGUOUS"), 
#     ct.c_int, ndpointer(ct.c_double, flags="C_CONTIGUOUS"), ndpointer(ct.c_double, flags="C_CONTIGUOUS"), 
#     ct.c_int, ndpointer(ct.c_int, flags="C_CONTIGUOUS")
# ]

def box_to_flexi_with_averaged_boundaries(xs, Xs, box, flx):
    shape = box.shape

    # ll ... lower left
    # tr ... top right
    lls, trs = flx.mesh.get_cell_coords()
    Is,Js,Ks = tuple(np.round((flx.Nout) * lls / flx.mesh.cellsize).astype(int).T)
    offsets  = ((Is * shape[1]) + Js) * shape[2] + Ks
    offsets  = np.require(offsets.ravel(), dtype=np.int32, requirements=['C', 'A'])

    xs  = np.require(xs.ravel(), dtype=np.double, requirements=['C', 'A'])
    Xs  = np.require(Xs.ravel(), dtype=np.double, requirements=['C', 'A'])
    box = np.require(box.ravel(), dtype=np.double, requirements=['C', 'A'])

    flxdata = np.empty(len(offsets) * len(Xs)**3, dtype=np.double)
    flxdata = np.require(flxdata.ravel(), dtype=np.double, requirements=['C', 'A', 'W'])

    lib.box_to_flexi_with_averaged_boundaries(
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
