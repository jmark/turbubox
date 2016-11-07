#!/usr/bin/env python3

import time
import flash
import hopr
import flexi
import gausslobatto
import scipy.interpolate
from scipy import ndimage
import numpy as np
import sys
import ulz
import interpolate

# just for testing
class FakeFlashFile:
    def __init__(self, domain, boxdims, fillby='constant', blockdims=None):
       	self.domain     = np.array(domain)
        self.gridsize   = np.array(boxdims)
        self.grid       = np.array([[0,0,0], self.gridsize-1])
        self.domainsize = np.abs(self.domain[1]-self.domain[0])
        self.cellsize   = self.domainsize / self.gridsize

        if blockdims:
            self.blocksize = np.array(blockdims)
        else:
            self.blocksize = None

        self.fillby = fillby
 
    @staticmethod
    def plateau(x,y,z):
        if np.abs(x) <= 0.2 and np.abs(y) <= 0.2 and np.abs(z) <= 0.2:
            return 1
        return 0

    @staticmethod
    def wiggle(X,Y,Z):
        return np.sin(4 * 2*np.pi * X) + np.sin(5 * 2*np.pi * Y) + np.sin(6 * 2*np.pi * Z)

    def data(self, dname):
        dom = self.domain.transpose()
        grd = self.gridsize
        ret = np.zeros(grd)
        fillby = self.fillby

        X,Y,Z = np.meshgrid(*tuple(ulz.mk_body_centered_linspace(*d, s) for d,s in zip(dom, grd)), indexing='ij')

        if dname == 'dens':

            if fillby == 'constant':
                ret[:] = 1.0

            elif fillby == 'planeX':
                ret[:] = X

            elif fillby == 'planeX+':
                ret = np.where(
                    (1/8-1/grd[0]<X) * (X<7/8+1/grd[0]) * \
                    (1/8-1/grd[1]<Y) * (Y<7/8+1/grd[1]) * \
                    (1/8-1/grd[2]<Z) * (Z<7/8+1/grd[2]), X,0*X)

            elif fillby == 'planeXYZ':
                ret = np.where(
                    (1/8-1/grd[0]<X) * (X<7/8+1/grd[0]) * \
                    (1/8-1/grd[1]<Y) * (Y<7/8+1/grd[1]) * \
                    (1/8-1/grd[2]<Z) * (Z<7/8+1/grd[2]), X+Y+Z,0*X)

            elif fillby == 'plane+wiggle':
                ret = np.where(
                    (1/16-1/grd[0]<X) * (X<15/16+1/grd[0]) * \
                    (1/16-1/grd[1]<Y) * (Y<15/16+1/grd[1]) * \
                    (1/16-1/grd[2]<Z) * (Z<15/16+1/grd[2]),
                    X+Y+Z + 0.5*self.wiggle(X+1/16,Y+1/16,Z+1/16),0*X)

            elif fillby == 'gaussianXYZ':
                ret = np.exp(-((X-0.5)**2 + (Y-0.5)**2 + (Z-0.5)**2)/2/0.02)

            elif fillby == 'stepsXYZ':
                ret = X + 100*X*(np.sin(4*2*np.pi*X) + 0.2*np.cos(20*2*np.pi*X) + 0.1*np.sin(20*2*np.pi*X))**2 

            else:
                raise NotImplementedError('unknow fillby: %s' % fillby)

        if dname == 'pres':
            # constant
                ret[:] = 1.0
            # plane
            #ret = X+Y

            # gaussian2d
            # ret = 1 + 1 * 10**1 * np.exp(-((X-0.5)**2 + (Y-0.5)**2)/2/0.02)

            # gaussian3d
            #return 1 + 1 * 10**1 * np.exp(-(X**2 + Y**2 + Z**2)/2/0.02)

            # sine3d
            #return np.abs(np.sin(10*(X**2 + Y**2 + Z**2)))

            # cube2d
            # return 1 + 2 * np.vectorize(plateau)(X,Y,0)

            # wiggle2d
            #ret = 1 + 10 * np.abs(np.sin((X-Y)/2/np.pi * 200) + np.sin((X+Y)/2/np.pi * 100))
        
            # sin2d
            #ret = 200 * (np.sin(X/np.pi/2) + np.cos(Y/np.pi/2))

            #ret = np.sin(2*np.pi*X)
   
            # sin2d
            #return X+Y+Z

        return ret
        
# =========================================================================== #

def linear_3d():
    sys.argv.reverse()
    progname = sys.argv.pop()
    flshfile = sys.argv.pop()
    meshfile = sys.argv.pop()
    flexfile = sys.argv.pop() 
    #method   = sys.argv.pop()

    if flshfile == '--generate':
        flash = FakeFlashFile([[0,0,0],[1,1,1]], [32]*3, [4]*3)
    else:
        flash = flash.File(flshfile)

    flexi = flexi.File(flexfile, hopr.CartesianMeshFile(meshfile))

    N = flexi.N

    nodes = mk_nodes_cuboid((N,N,N), 'gauss-lobatto')
    cells = flexi.mesh.get_cell_coords().transpose(1,0,2)
    rgrid = np.array([mk_cell_points(cell, nodes) for cell in cells]).reshape(-1, *(N+1,)*3, 3)
    cgrid = mk_body_centered_grid(flexi.mesh.domain, flash.gridsize, boundary=True)

    kappa = 5/3
    mu0   = 1
    prims = [] # primitive variables

    for dbname in 'dens velx vely velz pres magx magy magz'.split():
        prims.append(flash.data(dbname))

    cons    = [None]*len(prims)
    cons[0] = prims[0]            # density
    cons[1] = prims[0]*prims[1]   # momentum x
    cons[2] = prims[0]*prims[2]   # momentum y
    cons[3] = prims[0]*prims[3]   # momentum z
    cons[4] = prims[4]/(kappa-1) + prims[0]/2*(prims[1]**2+prims[2]**2+prims[3]**2) + (prims[5]**2+prims[6]**2+prims[7]**2)/2/mu0 # total energy
    cons[5] = prims[5]           # mag x
    cons[6] = prims[6]           # mag y
    cons[7] = prims[7]           # mag z

    method = 'linear'
    #method = 'nearest'
    #method = 'spline'
    order  = 1

    if method == 'spline':
        shape = rgrid.shape
        rgrid = rgrid.reshape(-1,3).T * (flash.gridsize[0]-1) - 1/flash.gridsize[0]/2
        #print(rgrid.shape)
        #print(rgrid)

    for nvar,convar in enumerate(cons):
        if method == 'linear' or method == 'nearest':
            burrito = wrap_in_guard_cells(convar)        
            itpl = scipy.interpolate.RegularGridInterpolator(cgrid, burrito, method=method)
        elif method == 'spline':
            itpl = lambda coords: ndimage.map_coordinates(convar, coords, order=order, mode='wrap').reshape(shape[:-1])
        else:
            raise NotImplemented("%s - unknown method" % method)

        idat = itpl(rgrid)
        print("Writing nvar '%d': (orig/itpl) min %f/%f max %f/%f" % \
                (nvar, convar.min(), idat.min(), convar.max(), idat.max()), file=sys.stderr)
        flexi.data[:,:,:,:,nvar] = idat.transpose(0,3,2,1)

def lagrange_3d_3rd_order():
    sys.argv.reverse()
    progname = sys.argv.pop()
    flshfile = sys.argv.pop()
    meshfile = sys.argv.pop()
    flexfile = sys.argv.pop() 

    flx = flexi.File(flexfile, hopr.CartesianMeshFile(meshfile))
    npoly = flx.npoly
    ntype = flx.nodetype
    gridsize = (flx.mesh.gridsize * (npoly+1)).astype(int)

    if '--generate=' in flshfile:
        fillby = flshfile.split('=')[-1]
        fls = FakeFlashFile([[0,0,0],[1,1,1]], gridsize, fillby=fillby)
    else:
        fls = flash.File(flshfile)

    xs   = ulz.mk_body_centered_linspace(-1,1,npoly+1)
    Xs   = ulz.mk_cartesian_product_3d(*[gausslobatto.mk_nodes(npoly, ntype)]*3)
    ipl  = gausslobatto.mk_lagrange_interpolator_3d(*[xs]*3,Xs)
    dens = fls.data('dens')
    pres = fls.data('pres')

    elemsize = flx.mesh.cellsize / (npoly+1)

    # ll = lower left
    # tr = top right
    for elemid,ll,tr in zip(range(0,flx.mesh.nrelems),*flx.mesh.get_cell_coords()):
        i,j,k    = tuple(np.round(ll/elemsize).astype(int))
        io,jo,ko = tuple(np.round(tr/elemsize).astype(int))

        flx.data[elemid,:,:,:,0] = ipl(dens[i:io,j:jo,k:ko]).reshape(*[npoly+1]*3).transpose(2,1,0)
        flx.data[elemid,:,:,:,1] = dens[i:io,j:jo,k:ko].transpose(2,1,0)
        flx.data[elemid,:,:,:,4] = ipl(pres[i:io,j:jo,k:ko]).reshape(*[npoly+1]*3).transpose(2,1,0)

    idat = flx.data[:,:,:,:,0]
    print("Writing nvar '%s': (orig/itpl) min %f/%f max %f/%f" % \
            ('dens', dens.min(), idat.min(), dens.max(), idat.max()), file=sys.stderr)

    idat = flx.data[:,:,:,:,4]
    print("Writing nvar '%s': (orig/itpl) min %f/%f max %f/%f" % \
            ('pres', pres.min(), idat.min(), pres.max(), idat.max()), file=sys.stderr)

def lagrange_3d_5th_order():
    sys.argv.reverse()
    progname = sys.argv.pop()
    flshfile = sys.argv.pop()
    meshfile = sys.argv.pop()
    flexfile = sys.argv.pop() 

    flx = flexi.File(flexfile, hopr.CartesianMeshFile(meshfile))
    npoly = flx.npoly
    ntype = flx.nodetype
    gridsize = (flx.mesh.gridsize * (npoly+1)).astype(int)

    if '--generate=' in flshfile:
        fillby = flshfile.split('=')[-1]
        fls = FakeFlashFile([[0,0,0],[1,1,1]], gridsize, fillby=fillby)
    else:
        fls = flash.File(flshfile)

    xs   = ulz.mk_body_centered_linspace(-1,1,npoly+1,withBoundaryNodes=True)
    #Xs   = ulz.mk_cartesian_product_3d(*[gausslobatto.mk_nodes(npoly, ntype)]*3)

    xs_  = ulz.mk_body_centered_linspace(-1,1,npoly+1)
    Xs   = ulz.mk_cartesian_product_3d(*[xs_]*3)
    ipl  = gausslobatto.mk_lagrange_interpolator_3d(*[xs]*3,Xs)

    dens    = fls.data('dens')
    burrito = ulz.wrap_in_guard_cells(dens)        
    elemsize = flx.mesh.cellsize / (npoly+1)

    start = time.time()
    # ll = lower left
    # tr = top right
    for elemid,ll,tr in zip(range(0,flx.mesh.nrelems),*flx.mesh.get_cell_coords()):
    #for elemid,ll,tr in zip(range(0,3),*flx.mesh.get_cell_coords()):
        i,j,k    = tuple(np.round(ll/elemsize).astype(int))
        io,jo,ko = tuple(np.round(tr/elemsize).astype(int))

        srcdata = burrito[i:io+2,j:jo+2,k:ko+2]
        snkdata = ipl(srcdata).reshape(*[npoly+1]*3)

        flx.data[elemid,:,:,:,0] = snkdata.transpose(2,1,0)

    print("Elapsed: %f s" % (time.time() - start))
    idat = flx.data[:,:,:,:,0]
    print("Writing nvar '%s': (orig/itpl) min %f/%f max %f/%f" % \
            ('dens', dens.min(), idat.min(), dens.max(), idat.max()), file=sys.stderr)


def lagrange_3d_5th_order2():
    sys.argv.reverse()
    progname = sys.argv.pop()
    flshfile = sys.argv.pop()
    meshfile = sys.argv.pop()
    flexfile = sys.argv.pop() 

    flx = flexi.File(flexfile, hopr.CartesianMeshFile(meshfile))
    npoly = flx.npoly
    ntype = flx.nodetype
    gridsize = (flx.mesh.gridsize * (npoly+1)).astype(int)

    if '--generate=' in flshfile:
        fillby = flshfile.split('=')[-1]
        fls = FakeFlashFile([[0,0,0],[1,1,1]], gridsize, fillby=fillby)
    else:
        fls = flash.File(flshfile)

    linspace = ulz.mk_body_centered_linspace(-1,1,npoly+1, withBoundaryNodes=True)
    xs,ys,zs = [linspace]*3

    Linspace = ulz.mk_body_centered_linspace(-1,1,npoly+1) 
    Xs,Ys,Zs = np.meshgrid(*[Linspace]*3, indexing='ij')

    ipl = interpolate.lagrange_interpolate_3d

    dens     = fls.data('dens')
    burrito  = ulz.wrap_in_guard_cells(dens)        
    elemsize = flx.mesh.cellsize / (npoly+1)

    # ll ... lower left
    # tr ... top right
    lls, trs = flx.mesh.get_cell_coords()
    snk = flx.data

    I,J,K    = tuple(np.round(lls/elemsize).astype(int).T)
    IO,JO,KO = tuple(np.round(trs/elemsize).astype(int).T + 2)
    shp      = (npoly+1,)*3

    #snkdata = flash_to_flexi(Is, IOs, xs, ys, zs, srcdata, Xs, Ys, Zs)
    start = time.time()

    for i,j,k,io,jo,ko,ll,tr,snk in zip(I,J,K,IO,JO,KO,lls,trs,snk):
        srcdata         = burrito[i:io,j:jo,k:ko]
        snkdata         = ipl(xs,ys,zs,srcdata,Xs,Ys,Zs).reshape(shp)
        snk[:,:,:,0]    = snkdata.transpose(2,1,0)

        #print("%d / %d" % (i,flx.mesh.nrelems))

        # print(srcdata.shape)
        # print(xs.shape)
        # print(ys.shape)
        # print(zs.shape)

    print("Elapsed: %f s" % (time.time() - start))
    idat = flx.data[:,:,:,:,0]
    print("Writing nvar '%s': (orig/itpl) min %f/%f max %f/%f" % \
            ('dens', dens.min(), idat.min(), dens.max(), idat.max()), file=sys.stderr)

def lagrange_3d_5th_order3():
    sys.argv.reverse()
    progname = sys.argv.pop()
    flshfile = sys.argv.pop()
    meshfile = sys.argv.pop()
    flexfile = sys.argv.pop() 

    flx = flexi.File(flexfile, hopr.CartesianMeshFile(meshfile))
    npoly = flx.npoly
    ntype = flx.nodetype
    gridsize = (flx.mesh.gridsize * (npoly+1)).astype(int)

    if '--generate=' in flshfile:
        fillby = flshfile.split('=')[-1]
        fls = FakeFlashFile([[0,0,0],[1,1,1]], gridsize, fillby=fillby)
    else:
        fls = flash.File(flshfile)

    xs  = ulz.mk_body_centered_linspace(-1,1,npoly+1, withBoundaryNodes=True)
    box = ulz.wrap_in_guard_cells(fls.data('dens'))        

    #xs  = ulz.mk_body_centered_linspace(-1,1,npoly+1)
    #box = fls.data('dens')

    #Xs  = ulz.mk_body_centered_linspace(-1,1,npoly+1) 
    Xs  = gausslobatto.mk_nodes(npoly, ntype)

    start = time.time()
    snkdata = interpolate.box_to_flexi(xs, Xs, box, flx)
    print("Elapsed: %f s" % (time.time() - start))

    flx.data[:,:,:,:,0] = snkdata.transpose(0,3,2,1)

    srcdata = box 
    idat = flx.data[:,:,:,:,0]
    print("Writing nvar '%s': (orig/itpl) min %f/%f max %f/%f" % \
            ('dens', srcdata.min(), idat.min(), srcdata.max(), idat.max()), file=sys.stderr)


if __name__ == '__main__':
    #linear_3d()
    #lagrange_3d_3rd_order()
    #lagrange_3d_5th_order()
    #lagrange_3d_5th_order2()
    lagrange_3d_5th_order3()
