import itertools
import hashlib
import numpy as np
import os

def gen_key(unit,grid,solv,deli='/'):
    return '%s%s%d%s%s' % (unit,deli,grid,deli,solv)

def gen_amr_path(key):
    return '/srv/data/FLASH/stirturb/mach-2.0/%s' % key

def gen_pickle_path(key):
    cachedir = '/tmp'
    #uniqueid = hashlib.md5(key.encode()).hexdigest()
    uniqueid = key.replace('/','.')
    return '%s/time-evolution.%s.pandas.pickle' % (cachedir,uniqueid)

def gen_hdf5_path(key):
    cachedir = '/tmp'
    uniqueid = key.replace('/','.')
    return '%s/time-evolution.%s.pandas.h5' % (cachedir,uniqueid)

units = 'unit cgs'.split()
grids = [16,24,32,48,64]
solvs = '8w b3 b5 es'.split()

def load_datafiles(units,grids,solvs):
    keys = [gen_key(*key) for key in itertools.product(units,grids,solvs)]
    dfs  = [pd.read_hdf(gen_hdf5_path(key),'table') for key in keys]

    return dict(zip(keys,dfs))

def turntime(key):
    if 'unit' in key:
        return 1/2.0
    elif 'cgs' in key:
        return 4.385e+14


def sort_unstructured_grid(ugrid, udata, doReshape=True):
    ndim = ugrid.shape[-1]
    tdat = ugrid.dtype

    # structured/sorted grid + sorting indices
    grid,idx = np.unique(ugrid.view(tdat.descr * ndim), return_index=True)

    grid = grid.view(tdat).reshape(-1,ndim)
    data = udata.ravel()[idx]

    if doReshape:
        # try to reshape sorted data to ndimensional block of n cells
        n = int(np.round(np.power(grid.shape[0],1/ndim)))     
        grid = grid.reshape([n]*ndim + [ndim])
        data = data.reshape([n]*ndim)

    return grid, data

def mk_cartesian_product_2d(xs,ys):
    return np.transpose(np.meshgrid(xs,ys)).reshape(-1,2)

def mk_cartesian_product_3d(xs,ys,zs):
    return np.roll(np.transpose(np.meshgrid(xs,ys,zs)).reshape(-1,3),1,axis=1)

def mk_body_centered_linspace(infimum, supremum, nNodes, withBoundaryNodes=False):
    """
    Make regular body centered linear space w/ or w/o neighboring boundary nodes.
    """

    domsize = np.abs(supremum - infimum)
    offset  = domsize / nNodes / 2

    if withBoundaryNodes:
        nNodes     = nNodes + 2
        infimum    = infimum  - offset
        supremum   = supremum + offset
    else:
        infimum    = infimum  + offset
        supremum   = supremum - offset

    return np.linspace(infimum, supremum, nNodes, endpoint=True)

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

def transform_to_ref_space(left, right, nodes):
    return left + (right-left) * (nodes+1)/2


def diff_x(f,Delta=1):
    return (np.roll(f,-1,axis=0) - np.roll(f,1,axis=0))/2./Delta

def diff_y(f,Delta=1):
    return (np.roll(f,-1,axis=1) - np.roll(f,1,axis=1))/2./Delta

def diff_z(f,Delta=1):
    return (np.roll(f,-1,axis=2) - np.roll(f,1,axis=2))/2./Delta

def curl(X,Y,Z,Dx,Dy,Dz):
    dX = (diff_y(Z,Dy) - diff_z(Y,Dz))
    dY = (diff_z(X,Dz) - diff_x(Z,Dx))
    dZ = (diff_x(Y,Dx) - diff_y(X,Dy))
    
    return (dX,dY,dZ)

def norm(X,Y,Z):
    return X**2 + Y**2 + Z**2

def find_file(fname, paths):
    for path in paths:
        for root, dirs, files in os.walk(path):
            if fname in files:
                return os.path.join(root, fname)
    raise FileNotFoundError("Cannot find '%s' in any of %s." % (fname, paths))
