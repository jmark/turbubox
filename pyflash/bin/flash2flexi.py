#!/usr/bin/env python3

import flash
import hopr
import flexi
import gausslobatto
import scipy.interpolate
import numpy as np
import sys

def mk_body_centered_grid(domain, shape, boundary=False):
    """
    Make regular body centered grid w/ periodic boundaries.
    
    Parameters
    ----------
    domain : ndarray of shape (2,ndim)
        lower left and upper right corners of the (cubiod) domain
    shape : ndarray of shape (ndim)
        grid/shape size
    
    Returns
    -------
    ndarray of shape (ndim, shape+2)
        n-dimensional linspace
    """
    domain = np.array(domain)
    shape  = np.array(shape)

    domainSize = np.abs(domain[1] - domain[0])
    offset = domainSize / shape / 2

    if boundary:
        newshape   = shape+2
        infimum    = domain[0] - offset
        supremum   = domain[1] + offset
    else:
        newshape = shape
        infimum    = domain[0] + offset
        supremum   = domain[1] - offset

    return np.array([np.linspace(inf, supr, shp, endpoint=True) for inf, supr, shp in zip(infimum, supremum, newshape)])

def wrap_in_guard_cells(stone):
    """
    Wrap a cubiod block of data with 'guard cells' which
    represents the opposite side - periodic boundaries.

    Parameters
    ----------
    stone : ndarray

    Returns
    -------
    ndarray of shape (stone[0]+2, stone[1]+2, ...)
    """
    plum = np.zeros(np.array(stone.shape)+2)
    
    # fill stone of plum
    plum[1:-1,1:-1,1:-1] = stone

    # wrap up stone with pulp 
    plum[ 0,1:-1,1:-1] = stone[-1,:,:]
    plum[-1,1:-1,1:-1] = stone[ 0,:,:]

    plum[1:-1, 0,1:-1] = stone[:,-1,:]
    plum[1:-1,-1,1:-1] = stone[:, 0,:]

    plum[1:-1,1:-1, 0] = stone[:,:,-1]
    plum[1:-1,1:-1,-1] = stone[:,:, 0]

    # don't forget the edges
    plum[:, 0, 0] = plum[:, 0, -2]
    plum[:,-1, 0] = plum[:,-1, -2]
    plum[:, 0,-1] = plum[:, 0,  1]
    plum[:,-1,-1] = plum[:,-1,  1]

    plum[ 0,:, 0] = plum[ 0, :,-2]
    plum[-1,:, 0] = plum[-1, :,-2]
    plum[ 0,:,-1] = plum[ 0, :, 1]
    plum[-1,:,-1] = plum[-1, :, 1]

    plum[ 0, 0,:] = plum[ 0,-2, :]
    plum[-1, 0,:] = plum[-1,-2, :]
    plum[ 0,-1,:] = plum[ 0, 1, :]
    plum[-1,-1,:] = plum[-1, 1, :]

    return plum

def mk_nodes_cuboid(Ns, ntype='gauss'): # shape: (N,N,N,3)
    """ Calculate a 3D array of gauss/gauss-lobatto nodes.

    Parameters
    ----------
    Ns : sequence of 3 polynomial degrees for x,y,z direction resp.
    ntype : string of demanded node type: 'gauss' / 'gauss-lobatto'

    Returns
    -------
    ndarray of shape (Nx*Ny*Nz, 3)
    """
    if ntype == 'gauss':
        fun = gausslobatto.LegendreGaussNodesAndWeights
    elif ntype == 'gauss-lobatto':
        fun = gausslobatto.LegendreGaussLobattoNodesAndWeights
    else:
        raise KeyError("unknown node type: '%s'" % ntype)

    xs,_ = fun(Ns[0])
    ys,_ = fun(Ns[1])
    zs,_ = fun(Ns[2])

    return np.array([[x,y,z] for x in xs for y in ys for z in zs])

def mk_cell_points(coords, nodes):
    return coords[0] + (coords[1]-coords[0]) * (nodes+1)/2

def prod(arr):
    res = 1
    for a in arr:
        res *= a
    return res

# just for testing
class FakeFlashStateFile:
    def __init__(self, domain, boxdims, blockdims):
       	self.domain     = np.array(domain)
        self.gridsize   = np.array(boxdims)
        self.blocksize  = np.array(blockdims)
        self.grid       = np.array([[0,0,0], self.gridsize-1])
        self.domainsize = np.abs(self.domain[1]-self.domain[0])
        self.cellsize   = self.domainsize / self.gridsize

    def box(self, dname):
        dom = self.domain.transpose()
        grd = self.gridsize

        #return np.arange(1,prod(self.gridsize)+1).reshape(self.gridsize)
        x,y,z = tuple(np.linspace(*d, g) for d,g in zip(dom,grd))
        X,Y,Z = np.meshgrid(x, y, z)

        return np.sin(X**2 + Y**2 + Z**2)

if __name__ == '__main__':
    sys.argv.reverse()
    progname = sys.argv.pop()
    flshfile = sys.argv.pop()
    meshfile = sys.argv.pop()
    flexfile = sys.argv.pop() 

    flash = flash.File(flshfile)
    flexi = flexi.File(flexfile, hopr.CartesianMeshFile(meshfile))

    N = flexi.N

    nodes = mk_nodes_cuboid((N,N,N), 'gauss-lobatto')
    cells = flexi.mesh.get_cell_coords().transpose(1,0,2)
    rgrid = np.array([mk_cell_points(cell, nodes) for cell in cells])
    rgrid = rgrid.reshape(rgrid.shape[0], N+1, N+1, N+1, 3)
    cgrid = mk_body_centered_grid(flexi.mesh.domain, flash.gridsize, boundary=True)

    for nvar, dbname in enumerate('dens velx vely velz pres magx magy magz'.split()):
        data = wrap_in_guard_cells(flash.data(dbname))
        itpl = scipy.interpolate.RegularGridInterpolator(cgrid, data, method='linear')
        flexi.data[:,:,:,:,nvar] = itpl(rgrid).transpose(0,3,2,1)
        print("'%s' successfully written!" % dbname, file=sys.stderr)
