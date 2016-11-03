import numpy as np
import ctypes as ct
from numpy.ctypeslib import ndpointer

libpath = '/home/jmark/projects/stirturb/turbubox/pyflash/lib/libinterpolate.so'
lib     = ct.cdll.LoadLibrary(libpath)

lib.lagrange_interpolate_3d.argtypes = [
    # xs,ys,zs,fs
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"), ct.c_int,
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"), ct.c_int,
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"), ct.c_int,
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"), ct.c_int,

    # Xs,Ys,Zs,Fs
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"), ct.c_int]

lib.flash_to_flexi.argtypes = [
    # I,J,K
    ndpointer(ct.c_int, flags="C_CONTIGUOUS"),
    ndpointer(ct.c_int, flags="C_CONTIGUOUS"),
    ndpointer(ct.c_int, flags="C_CONTIGUOUS"),

    # IO,JO,KO
    ndpointer(ct.c_int, flags="C_CONTIGUOUS"),
    ndpointer(ct.c_int, flags="C_CONTIGUOUS"),
    ndpointer(ct.c_int, flags="C_CONTIGUOUS"), ct.c_int,

    # xs,ys,zs,fss
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"), ct.c_int,
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"), ct.c_int,
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"), ct.c_int,
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"), ct.c_int, ct.c_int, ct.c_int, ct.c_int,

    # Xs,Ys,Zs,Fss
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"), ct.c_int]

def lagrange_interpolate_3d(xs,ys,zs,fs,Xs,Ys,Zs):
    xs = np.require(xs.ravel(), dtype=np.double, requirements=['C', 'A'])
    ys = np.require(ys.ravel(), dtype=np.double, requirements=['C', 'A'])
    zs = np.require(zs.ravel(), dtype=np.double, requirements=['C', 'A'])
    fs = np.require(fs.ravel(), dtype=np.double, requirements=['C', 'A'])

    xslen = len(xs)
    yslen = len(ys)
    zslen = len(zs)
    fslen = len(fs)

    assert xslen * yslen * zslen == fslen

    Xs = np.require(Xs.ravel(), dtype=np.double, requirements=['C', 'A'])
    Ys = np.require(Ys.ravel(), dtype=np.double, requirements=['C', 'A'])
    Zs = np.require(Zs.ravel(), dtype=np.double, requirements=['C', 'A'])

    Fs = np.empty_like(Xs)
    Fs = np.require(Fs.ravel(), dtype=np.double, requirements=['C', 'A', 'W'])

    Fslen = len(Fs)
    assert len(Xs) == Fslen
    assert len(Ys) == Fslen
    assert len(Zs) == Fslen

    lib.lagrange_interpolate_3d(
        xs,xslen,
        ys,yslen,
        zs,zslen,
        fs,fslen,
        Xs,Ys,Zs,Fs,Fslen)

    return Fs 


def flash_to_flexi(Is,Js,Ks, IOs,JOs,KOs, xs,ys,zs,fss, Xs,Ys,Zs):

    Is = np.require(Is.ravel(), dtype=np.int32, requirements=['C', 'A'])
    Js = np.require(Js.ravel(), dtype=np.int32, requirements=['C', 'A'])
    Ks = np.require(Ks.ravel(), dtype=np.int32, requirements=['C', 'A'])

    IOs = np.require(IOs.ravel(), dtype=np.int32, requirements=['C', 'A'])
    JOs = np.require(JOs.ravel(), dtype=np.int32, requirements=['C', 'A'])
    KOs = np.require(KOs.ravel(), dtype=np.int32, requirements=['C', 'A'])

    xs = np.require(xs.ravel(), dtype=np.double, requirements=['C', 'A'])
    ys = np.require(ys.ravel(), dtype=np.double, requirements=['C', 'A'])
    zs = np.require(zs.ravel(), dtype=np.double, requirements=['C', 'A'])

    Xs = np.require(Xs.ravel(), dtype=np.double, requirements=['C', 'A'])
    Ys = np.require(Ys.ravel(), dtype=np.double, requirements=['C', 'A'])
    Zs = np.require(Zs.ravel(), dtype=np.double, requirements=['C', 'A'])

    fssshape = fss.shape
    fss = np.require(fss.ravel(), dtype=np.double, requirements=['C', 'A'])

    IJKslen = len(Is)
    xslen = len(xs)
    yslen = len(ys)
    zslen = len(zs)
    fslen = xslen*yslen*zslen
    Fslen = len(Xs)

    Fss = np.empty(IJKslen * Fslen)
    Fss = np.require(Fss.ravel(), dtype=np.double, requirements=['C', 'A', 'W'])

    lib.flash_to_flexi(
        Is,Js,Ks, IOs,JOs,KOs, IJKslen, 
        xs,xslen,
        ys,yslen,
        zs,zslen,
        fss,*fssshape,fslen,
        Xs,Ys,Zs,Fss,Fslen)

    return Fss
