import numpy as np
import ctypes as ct
from numpy.ctypeslib import ndpointer
import sys
import os

def find_file(fname, paths):
    for path in paths:
        for root, dirs, files in os.walk(path):
            if fname in files:
                return os.path.join(root, fname)
    raise FileNotFoundError("Cannot find '%s' in any of %s." % (fname, paths))

lib = ct.cdll.LoadLibrary(find_file('libfortinterpolate.so', sys.path))

ptr_int8    = ndpointer(ct.c_int8,      flags="C_CONTIGUOUS")
ptr_int32   = ndpointer(ct.c_int32,     flags="C_CONTIGUOUS")
ptr_double  = ndpointer(ct.c_double,    flags="C_CONTIGUOUS")

def carray(ndarray, dtype=None):
    return np.require(ndarray, dtype=dtype, requirements=['C','A'])

lib.foo.argtypes = (
    ct.c_int32, ct.c_int32, ptr_double, ptr_double,
)

def foo(input):
    output = np.zeros_like(input)

    lib.foo(
        input.shape[0], input.shape[1], carray(input), carray(output),
    )

    return output

if False:

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
    t_nint    = ndpointer(ct.c_int, flags="C_CONTIGUOUS")
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

        lib.box_to_elements(box.shape[0], box.shape[1], box.shape[2], boxptr, flx.mesh.nrelems, ofsptr, elems[0].shape[0], elems[0].shape[1], elems[0].shape[2], elemsptr, neighbors)

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

        lib.elements_to_box(box.shape[0], box.shape[1], box.shape[2], boxptr, mesh.nrelems, ofsptr, elems[0].shape[0], elems[0].shape[1], elems[0].shape[2], elemsptr)

        return boxptr.reshape(box.shape)

    # =========================================================================== #

    lib.elements_to_box_fv.argtypes = [
        t_int, t_int, t_int, t_ndouble,
        t_int, t_ndouble, 
        t_int, t_int, t_int, t_ndouble,
        t_nint
    ]

    def elements_to_box_fv(elems, mesh, box, fvs):
        # lower left corners normed to unit intervall
        lls      = (mesh.domain[0] + mesh.elemcoords[0])/mesh.domainsize

        #box      = np.zeros(elems[0].shape * mesh.meshshape)
        offsets  = np.array(box.shape) * lls

        boxptr   = np.require(box.ravel(), dtype=np.double, requirements=['C','A'])
        elemsptr = np.require(elems.ravel(), dtype=np.double, requirements=['C','A'])
        ofsptr   = np.require(offsets.ravel(), dtype=np.double, requirements=['C', 'A'])
        fvsptr   = np.require(fvs.ravel(), dtype=np.int, requirements=['C','A'])

        lib.elements_to_box_fv(
            box.shape[0], box.shape[1], box.shape[2], boxptr,
            mesh.nrelems, ofsptr,
            elems[0].shape[0], elems[0].shape[1], elems[0].shape[2], elemsptr,
            fvs)

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

        lib.box_to_elements_avg_boundaries(box.shape[0], box.shape[1], box.shape[2], boxptr, flx.mesh.nrelems, ofsptr, elems[0].shape[0], elems[0].shape[1], elems[0].shape[2], elemsptr)

        return elemsptr.reshape(elems.shape)

    # =========================================================================== #
    # void
    # change_basis_3d(
    #     const int nelems, const int nn,
    #     const double *Vdm, const double *fss, double *Fss);

    lib.change_basis_3d.argtypes = [
       t_int, t_int, t_ndouble, t_ndouble, t_ndouble
    ]

    def change_basis(Vd,fs):
        nn,NN = Vd.shape
        Fs = np.empty([len(fs)]+3*[NN])

        Vdptr = np.require(Vd, dtype=np.double, requirements=['C','A'])
        fsptr = np.require(fs, dtype=np.double, requirements=['C','A'])
        Fsptr = np.require(Fs, dtype=np.double, requirements=['C','A','W'])

        lib.change_basis_3d(
            len(fs), NN, Vdptr, fsptr, Fsptr
        )

        return Fsptr.reshape(Fs.shape)

    # =========================================================================== #
    # void
    # change_basis_3d_2(
    #     const int nelems, const int nn,
    #     const double *Vdm, const double *fss, double *Fss);

    lib.change_basis_3d_2.argtypes = [
       t_int, t_int, t_ndouble, t_ndouble, t_ndouble
    ]

    def change_basis_2(Vd,fs):
        nn,NN = Vd.shape
        Fs = np.empty([len(fs)]+3*[NN])

        Vdptr = np.require(Vd, dtype=np.double, requirements=['C','A'])
        fsptr = np.require(fs, dtype=np.double, requirements=['C','A'])
        Fsptr = np.require(Fs, dtype=np.double, requirements=['C','A','W'])

        lib.change_basis_3d_2(
            len(fs), NN, Vdptr, fsptr, Fsptr
        )

        return Fsptr.reshape(Fs.shape)

    # =========================================================================== #
    # void
    # change_grid_space_2d_2(
    #     const int nelems,
    #     const int nx, const int ny,
    #     const int Nx, const int Ny,
    #     const double *Lss, const double *fss, double *Fss);

    lib.change_grid_space_2d_2.argtypes = [
        t_int,
        t_int, t_int,
        t_int, t_int,
        t_ndouble, t_ndouble, t_ndouble
    ]

    def change_grid_space_2d_2(Ls,fs):
        sh = Ls.shape
        Fs = np.empty([len(fs), sh[0], sh[1]])

        Lsptr = np.require(Ls, dtype=np.double, requirements=['C','A'])
        fsptr = np.require(fs, dtype=np.double, requirements=['C','A'])
        Fsptr = np.require(Fs, dtype=np.double, requirements=['C','A','W'])

        lib.change_grid_space_2d_2(
            len(fs), sh[2], sh[3], sh[0], sh[1], Lsptr, fsptr, Fsptr
        )

        return Fsptr.reshape(Fs.shape)

    # =========================================================================== #
    # void
    # change_grid_space_2d(
    #     const int nelems,
    #     const int nx, const int ny, const double *xs, double *fss, 
    #     const int Nx, const int Ny, const double *Xs, double *Fss);

    lib.change_grid_space_2d.argtypes = [
       t_int,
       t_int, t_int, t_ndouble, t_ndouble,
       t_int, t_int, t_ndouble, t_ndouble
    ]

    def change_grid_space_2d(fs,xs,Xs):
        Fs = np.empty([len(fs), len(Xs), len(Xs)])

        xsptr = np.require(xs, dtype=np.double, requirements=['C','A'])
        Xsptr = np.require(Xs, dtype=np.double, requirements=['C','A'])
        fsptr = np.require(fs, dtype=np.double, requirements=['C','A'])
        Fsptr = np.require(Fs, dtype=np.double, requirements=['C','A','W'])

        lib.change_grid_space_2d(
            len(fs),
            len(xs), len(xs), xsptr, fsptr,
            len(Xs), len(Xs), Xsptr, Fsptr,
        )

        return Fsptr

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
    # void
    # change_grid_space_dg_fv(
    #     const int nelems,
    #     const int nx, const int ny, const int nz ,const double *xs, double *fss, 
    #     const int Nx, const int Ny, const int Nz ,const double *Xs, double *Fss,
    #     const int *fvs);

    lib.change_grid_space_dg_fv.argtypes = [
        t_int,
        t_int, t_int, t_int, t_ndouble, t_ndouble,
        t_int, t_int, t_int, t_ndouble, t_ndouble,
        t_nint
    ]

    def change_grid_space_dg_fv(fs,xs,Xs,FV):
        Fs = np.empty([len(fs), len(Xs), len(Xs), len(Xs)])

        xsptr = np.require(xs, dtype=np.double, requirements=['C','A'])
        Xsptr = np.require(Xs, dtype=np.double, requirements=['C','A'])
        fsptr = np.require(fs, dtype=np.double, requirements=['C','A'])
        Fsptr = np.require(Fs, dtype=np.double, requirements=['C','A','W'])
        FVptr = np.require(FV, dtype=np.int32,  requirements=['C','A'])

        lib.change_grid_space_dg_fv(
            len(fs),
            len(xs), len(xs), len(xs), xsptr, fsptr,
            len(Xs), len(Xs), len(Xs), Xsptr, Fsptr,
            FVptr
        )

        return Fsptr.reshape(Fs.shape)

    # =========================================================================== #
    # void
    # change_grid_space_fv_dg(
    #     const int nelems,
    #     const int nx, const int ny, const int nz ,const double *xs, double *fss, 
    #     const int Nx, const int Ny, const int Nz ,const double *Xs, double *Fss,
    #     const int *fvs);

    lib.change_grid_space_fv_dg.argtypes = [
        t_int,
        t_int, t_int, t_int, t_ndouble, t_ndouble,
        t_int, t_int, t_int, t_ndouble, t_ndouble,
        t_nint
    ]

    def change_grid_space_fv_dg(fs,xs,Xs,FV):
        Fs = np.empty([len(fs), len(Xs), len(Xs), len(Xs)])

        xsptr = np.require(xs, dtype=np.double, requirements=['C','A'])
        Xsptr = np.require(Xs, dtype=np.double, requirements=['C','A'])
        fsptr = np.require(fs, dtype=np.double, requirements=['C','A'])
        Fsptr = np.require(Fs, dtype=np.double, requirements=['C','A','W'])
        FVptr = np.require(FV, dtype=np.int32,  requirements=['C','A'])

        lib.change_grid_space_fv_dg(
            len(fs),
            len(xs), len(xs), len(xs), xsptr, fsptr,
            len(Xs), len(Xs), len(Xs), Xsptr, Fsptr,
            FVptr
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
            shape[0], shape[1], shape[2],  box,
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
            shape[0], shape[1], shape[2], box,
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
            shape[0], shape[1], shape[2], box)

        return box.reshape(shape)

    # =========================================================================== #

    # void
    # blocks_to_box(
    #     const int  rlevel,  const int nblocks,
    #     const int *rlevels, const double *coords, const double *domain,
    #     const int nx, const int ny, const int nz, const double *blocks,
    #     const int Nx, const int Ny, const int Nz, const double *box,
    # );

    lib.blocks_to_box.argtypes = [
        ct.c_int, ct.c_int, 
        ndpointer(ct.c_int,    flags="C_CONTIGUOUS"),
        ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
        ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
        ct.c_int, ct.c_int, ct.c_int, ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
        ct.c_int, ct.c_int, ct.c_int, ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ]

    def blocks_to_box(rlevel, rlevels, coords, domain, blocks, box):
        _rlevels = np.require(np.ravel(rlevels), dtype=np.int32,  requirements=['C', 'A'])
        _coords  = np.require(np.ravel( coords), dtype=np.double, requirements=['C', 'A'])
        _domain  = np.require(np.ravel( domain), dtype=np.double, requirements=['C', 'A'])
        _blocks  = np.require(np.ravel( blocks), dtype=np.double, requirements=['C', 'A'])
        _box     = np.require(np.ravel(    box), dtype=np.double, requirements=['C', 'A'])

        lib.blocks_to_box(
            rlevel, len(_rlevels), _rlevels, _coords, _domain,
            blocks.shape[1], blocks.shape[2], blocks.shape[3], _blocks,
               box.shape[0],    box.shape[1],    box.shape[2], _box)

        return _box.reshape(box.shape)

    # =========================================================================== #
    # =========================================================================== #
    # =========================================================================== #

    ptr_int8    = ndpointer(ct.c_int8,      flags="C_CONTIGUOUS")
    ptr_int32   = ndpointer(ct.c_int32,     flags="C_CONTIGUOUS")
    ptr_double  = ndpointer(ct.c_double,    flags="C_CONTIGUOUS")

    def carray(ndarray, dtype=None):
        return np.require(ndarray, dtype=dtype, requirements=['C','A'])

    lib.morton_to_coords.argtypes = (
        ptr_int32, ptr_int8, 
        ptr_int32, ptr_int32,
        ptr_int32, ptr_double,
    )

    def morton_to_coords(levels, morton):
        coords = np.zeros((levels.shape[0],3))

        lib.morton_to_coords(
            carray(levels.shape, dtype=np.int32), carray(levels),
            carray(morton.shape, dtype=np.int32), carray(morton),
            carray(coords.shape, dtype=np.int32), carray(coords),
        )

        return coords

    lib.cells_to_image.argtypes = (
        ptr_int32, ptr_int8, 
        ptr_int32, ptr_int32,
        ptr_int32, ptr_double,
        ptr_int32, ptr_double,
        ct.c_int32
    )

    def cells_to_image(levels, morton, cells, image, method='nearest', gridlines=0):
        methods = dict(
            nearest  = 0,
            bilinear = 1,
            bicosine = 2,
        )

        if len(cells.shape) < 2:
            cshape = (cells.shape[0], 1,1)
     
            lib.cells_to_image(
                carray(levels.shape, dtype=np.int32), carray(levels),
                carray(morton.shape, dtype=np.int32), carray(morton),
                carray(      cshape, dtype=np.int32), carray(cells),
                carray( image.shape, dtype=np.int32), carray(image),
                methods[str.lower(method)], gridlines
            )
        else:
            lib.cells_to_image(
                carray(levels.shape, dtype=np.int32), carray(levels),
                carray(morton.shape, dtype=np.int32), carray(morton),
                carray( cells.shape, dtype=np.int32), carray(cells),
                carray( image.shape, dtype=np.int32), carray(image),
                methods[str.lower(method)], gridlines
            )

    lib.cells_to_image_3d.argtypes = (
        ptr_int32, ptr_int8, 
        ptr_int32, ptr_int32,
        ptr_int32, ptr_double,
        ptr_int32, ptr_double,
    )

    def cells_to_image_3d(levels, morton, cells, image):
        if len(cells.shape) < 2:
            cshape = (cells.shape[0], 1,1,1)
            lib.cells_to_image_3d(
                carray(levels.shape, dtype=np.int32), carray(levels),
                carray(morton.shape, dtype=np.int32), carray(morton),
                carray(      cshape, dtype=np.int32), carray(cells),
                carray( image.shape, dtype=np.int32), carray(image)
            )
        else:
            lib.cells_to_image_3d(
                carray(levels.shape, dtype=np.int32), carray(levels),
                carray(morton.shape, dtype=np.int32), carray(morton),
                carray( cells.shape, dtype=np.int32), carray(cells),
                carray( image.shape, dtype=np.int32), carray(image)
            )


    lib.cells_to_image.argtypes = (
        ptr_int32, ptr_int8, 
        ptr_int32, ptr_int32,
        ptr_int32, ptr_double,
        ptr_int32, ptr_double,
        ct.c_int32
    )

    # =========================================================================== #

    lib.cells_to_image_flash_ug_2d.argtypes = (
        ptr_int32, ptr_double,
        ptr_int32, ptr_double,
        ptr_int32, ptr_double,
        ptr_int32, ptr_double,
        ct.c_int32,
    )

    def cells_to_image_flash_ug_2d(coords, bsizes, blocks, image, method='nearest'):
        methods = dict(
            nearest  = 0,
            bilinear = 1,
            bicosine = 2,
        )

        lib.cells_to_image_flash_ug_2d(
            carray(coords.shape, dtype=np.int32), carray(coords),
            carray(bsizes.shape, dtype=np.int32), carray(bsizes),
            carray(blocks.shape, dtype=np.int32), carray(blocks),
            carray( image.shape, dtype=np.int32), carray(image),
            methods[str.lower(method)]
        )

    # =========================================================================== #

    lib.cells_to_image_titanic_patch_2d.argtypes = (
        ptr_int32, ptr_double,
        ptr_int32, ptr_double,
        ct.c_int32,
    )

    def cells_to_image_titanic_patch_2d(blocks, image, method='nearest'):
        methods = dict(
            nearest  = 0,
            bilinear = 1,
        )

        lib.cells_to_image_titanic_patch_2d(
            carray(blocks.shape, dtype=np.int32), carray(blocks),
            carray( image.shape, dtype=np.int32), carray(image),
            methods[str.lower(method)]
        )
