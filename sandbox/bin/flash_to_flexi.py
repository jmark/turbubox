#!/usr/bin/env pyturbubox

# This script is a mess and needs clean up. But it kinda works, so...

import time
import flash
import hopr
import flexi
import gausslobatto
import scipy.interpolate
from scipy import ndimage
import periodicbox
import numpy as np
import sys
import ulz
import interpolate
import dslopts
import pathlib
import scipy.ndimage
               
# =========================================================================== #

def linear_3d():
    sys.argv.reverse()
    progname = sys.argv.pop()
    flshfile = sys.argv.pop()
    meshfile = sys.argv.pop()
    flexfile = sys.argv.pop() 
    #method   = sys.argv.pop()

    if flshfile == '--generate':
        flash = flash.FakeFile([[0,0,0],[1,1,1]], [32]*3, [4]*3)
    else:
        flash = flash.File(flshfile)

    flexi = flexi.File(flexfile, hopr.CartesianMeshFile(meshfile))

    N = flexi.N

    nodes = mk_nodes_cuboid((N,N,N), 'gauss-lobatto')
    cells = flexi.mesh.get_cell_coords().transpose(1,0,2)
    rgrid = np.array([mk_cell_points(cell, nodes) for cell in cells]).reshape(-1, N+1, N+1, N+1, 3)
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
        fls = flash.FakeFile([[0,0,0],[1,1,1]], gridsize, fillby=fillby)
    else:
        fls = flash.File(flshfile)

    xs   = ulz.mk_body_centered_linspace(-1,1,npoly+1)
    Xs   = ulz.mk_cartesian_product_3d(*[gausslobatto.mk_nodes(npoly, ntype)]*3)
    ipl  = gausslobatto.mk_lagrange_interpolator_3d(xs,xs,xs,Xs)
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
        fls = flash.FakeFile([[0,0,0],[1,1,1]], gridsize, fillby=fillby)
    else:
        fls = flash.File(flshfile)

    xs   = ulz.mk_body_centered_linspace(-1,1,npoly+1,withBoundaryNodes=True)
    #Xs   = ulz.mk_cartesian_product_3d(*[gausslobatto.mk_nodes(npoly, ntype)]*3)

    xs_  = ulz.mk_body_centered_linspace(-1,1,npoly+1)
    Xs   = ulz.mk_cartesian_product_3d(*[xs_]*3)
    ipl  = gausslobatto.mk_lagrange_interpolator_3d(xs,xs,xs,Xs)

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

methods = []

def method(box, flx, dname=None):
    "Without neighboring cells: n-th order interpolation (old version)"
    xs  = ulz.mk_body_centered_linspace(-1,1, flx.Nout)
    Xs  = gausslobatto.mk_nodes(flx.Nout-1, flx.nodetype) # target grid space

    return interpolate.box_to_flexi(xs, Xs, box, flx)

methods.append(method)

def method(box, flx, dname=None):
    "With neighboring cells: (n+2)-th order interpolation (old version)"
    xs  = ulz.mk_body_centered_linspace(-1,1,flx.Nout, withBoundaryNodes=True)
    Xs  = gausslobatto.mk_nodes(flx.Nout-1, flx.nodetype) # target grid space

    return interpolate.box_to_flexi(xs, Xs, ulz.wrap_in_guard_cells(box), flx)

methods.append(method)

def method(box, flx, dname=None):
    "Like '0' (new version)"
    els = interpolate.box_to_elements(box,flx)

    xs  = ulz.mk_body_centered_linspace(-1,1, els.shape[1])
    Xs  = gausslobatto.mk_nodes(flx.Nout-1, flx.nodetype) # target grid space

    return interpolate.change_grid_space(els, xs, Xs)

methods.append(method)

def method(box, flx, dname=None):
    "Like '1' (new version)"
    els = interpolate.box_to_elements(box, flx, neighbors = 1)

    xs  = ulz.mk_body_centered_linspace(-1,1, els.shape[1] - 2, withBoundaryNodes=True)
    Xs  = gausslobatto.mk_nodes(flx.Nout-1, flx.nodetype) # target grid space

    return interpolate.change_grid_space(els, xs, Xs)

methods.append(method)

def method(box, flx, dname=None):
    "like previous one + various filters"

    import scipy.ndimage.filters
    import scipy.ndimage.interpolation

    def pipe(source, *filters):
        sink = source
        for filter in filters:
            sink = filter(sink) 
        return sink

    def blur(src):
        factor = 2
        return scipy.ndimage.filters.gaussian_filter(src, factor)

    def zoom(src):
        factor = 1
        #factor = 1/2
        #factor = 396 / 256
        return scipy.ndimage.interpolation.zoom(src, factor, order=1, mode='wrap')

    def scale(src):
        factor = 5
        return factor * src

    if dname in 'dens pres':
        avg = box.mean()
        box = pipe(box,blur,zoom)
        box = avg / np.mean(box) * box

    if dname in 'velx vely velz':
        box = pipe(box,blur,scale,zoom)

    return methods[3](box, flx)

methods.append(method)

def method(box, flx, dname=None):
    "With neighboring cells which get averaged with boundary cells: n-th order interpolation."
    els = interpolate.box_to_elements_avg_boundaries(box, flx)

    xs  = ulz.mk_body_centered_linspace(-1,1, els.shape[1])
    Xs  = gausslobatto.mk_nodes(flx.Nout-1, flx.nodetype)
    xs[0] = -1; xs[-1]  =  1

    return interpolate.change_grid_space(els, xs, Xs)
        
methods.append(method)

def method(box, flx, dname=None):
    "Just copy over data. No interpolation!"
    return interpolate.box_to_elements(box, flx, neighbors = 0)
       
methods.append(method)

def lagrange_3d_5th_order3():
    def ExistingPath(arg):
        if '--generate=' in arg:
            fillby = arg.split('=',1)[-1]
            return fillby

        pth = pathlib.Path(arg)
        if not pth.exists():
            raise OSError("'%s' does not exists!" % pth)
        return pth

    def ConstrainedInt(arg):
        nr = int(arg)
        if 0 <= nr <= len(methods)-1:
            return nr
        raise ValueError("Method number must be within 0 and 4: '%d' given!")

    appendix = "  Methods are:\n\n" + "\n".join("    " + str(i) + " -> " + x.__doc__ for i,x in enumerate(methods))

    with dslopts.Manager(scope=globals(),appendix=appendix) as mgr:
        mgr.add(name='flashfile' ,desc='flash file path' ,type=ExistingPath)
        #mgr.add(name='meshfile'  ,desc=' mesh file path' ,type=ExistingPath)
        mgr.add(name='flexifile' ,desc='flexi file path' ,type=ExistingPath)
        mgr.add(name='method'    ,desc='method nr: 0-4'  ,type=ConstrainedInt, default=3)

    flx = periodicbox.File(str(flexifile), mode='r+')

    if isinstance(flashfile, pathlib.Path):
        fls = flash.File(str(flashfile))
    else:
        fls = flash.FakeFile([[0,0,0],[1,1,1]],flx.mesh.gridsize * flx.Nout, fillby=flashfile)

    # interpolate
    prims = []
    print(" Primitive vars:")
    print("")
    print("  var   |       min    ->    min     |       max    ->    max     ")
    print("  ------|----------------------------|----------------------------")

    #for dbname in 'dens velx vely velz pres magx magy magz'.split():
    for dbname in 'dens velx vely velz pres'.split():
        box = fls.data(dbname)
        els = methods[method](box, flx, dbname)

        print("  %s  | % 12.5f % 12.5f  | % 12.5f % 12.5f" % \
                (dbname, box.min(), els.min(), box.max(), els.max()), file=sys.stderr)

        prims.append(els)

    print("")
    print(" Conservative vars:")
    print("")
    print("  var   |       min     |     max     ")
    print("  ------|---------------|-------------")

    cons = ulz.navier_primitive_to_conservative(prims)
    # write to file
    for i,con in enumerate(cons):
        flx.data[:,:,:,:,i] = con.transpose(0,3,2,1)
        print("   % 2d   | % 12.5f  | % 12.5f" % (i, flx.data[:,:,:,:,i].min(), flx.data[:,:,:,:,i].max()), file=sys.stderr)

    print("")
    print(" Primitive vars:")
    print("")
    print("  var   |       min     |     max     ")
    print("  ------|---------------|-------------")

    for i,prim in enumerate(ulz.navier_conservative_to_primitive(cons)):
        print("   % 2d   | % 12.5f  | % 12.5f" % (i, prim.min(), prim.max()), file=sys.stderr)

if __name__ == '__main__':
    #linear_3d()
    #lagrange_3d_3rd_order()
    #lagrange_3d_5th_order()
    #lagrange_3d_5th_order2()
    lagrange_3d_5th_order3()
