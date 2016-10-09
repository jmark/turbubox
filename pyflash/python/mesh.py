"""
FLASH Code, FLEXI, HOPR utilities
"""

import h5py as h5
import numpy as np
import sys
import scipy.interpolate
import interpol

def dataset_to_dict(dataset):
    """
    Transform a key/value dataset to a pythonic dictionary

    Parameters
    ----------
    dataset : str
        dataset name

    Returns
    -------
    dictionary
    """
    return dict([(k.strip().decode(),v) for (k,v) in dataset])

def rbcgwpb(domain, shape):
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
    ndarray of shape (ndim, xx)
        n-dimensional linspace
    """

    domainSize = np.abs(domain[1] - domain[0])
    cellsize   = domainSize / shape
    infimum    = domain[0] - cellsize/2
    supremum   = domain[1] + cellsize/2

    return np.array([np.linspace(infimum[i], supremum[i], shape[i], endpoint=True) for i in range(0,3)])

def add_periodic_boundaries(stone):
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

def node_cube(N, ntype='gauss'): # shape: (N,N,N,3)
    if ntype == 'gaus':
        xs,ws = interpol.LegendreGaussNodesAndWeights(N)
    elif ntype == 'gauss-lobatto':
        xs,ws = interpol.LegendreGaussLobattoNodesAndWeights(N)
    else:
        raise KeyError("unknown node type: '%s'" % ntype)

    return np.array([[x,y,z] for x in xs for y in xs for z in xs])

def nodes_to_points(cell, nodes):
    return cell[0] + (cell[0] + (cell[1]-cell[0])/2) * nodes

def flash2flexi(flash, flexi):
    nodes = node_cube(flexi.pN, flexi.pT)

    for nvar, dbname in enumerate('density velx vely velz press bx by bz'.split()):
        box  = add_periodic_boundaries(flash.get(dbname))
        grid = rbcgwpb(flexi.hopr.domain, np.array(box.shape))
        itpl = scipy.interpolate.RegularGridInterpolator(grid, box, method='linear')
        
        for cellid, cell in flexi.hopr.cells:
            points = nodes_to_points(cell, nodes) # shape: (N,N,N,3)
            flexi.data[cellid,:,:,:,nvar] = itpl(points)

def prod(arr):
    res = 1
    for a in arr:
        res *= a
    return res

class FakeAMR:

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

        return np.arange(1,prod(self.gridsize)+1).reshape(self.gridsize)
        x,y,z = tuple(np.linspace(*d, g) for d,g in zip(dom,grd))
        X,Y,Z = np.meshgrid(x, y, z)

        return np.sin(X**2 + Y**2 + Z**2)

## --- Classes --- ##

class H5File(h5.File):

    def __init__(self,fpath,mode):
        try:
            super().__init__(fpath,mode)
        except OSError as e:
            print(sys.stderr,"OSError: ", e)
            sys.exit(1)
    
    def get(self,dname):
        db = super(H5File,self).get(dname)
        if db:
            return db
        else:
            raise KeyError(dname)

class AMR:

    def __init__(self,fpath):
        self.h5file = H5File(fpath,'r')

        self.simInfo        = self.get('sim info')

        self.integerScalars = dataset_to_dict(self.get('integer scalars'))
        self.realScalars    = dataset_to_dict(self.get('real scalars'))
        self.integerRuntime = dataset_to_dict(self.get('integer runtime parameters'))
        self.realRuntime    = dataset_to_dict(self.get('real runtime parameters'))

        self.gridSize = np.array([
            self.integerScalars[x]*2**(self.maxRefineLevel-1) for x in 'nxb nyb nzb'.split()])

        if str(self.simInfo['setup call'][0]).find('+ug') >= 0:
            self.gridSize *= 2 

        self.grid = np.array([0,0,0], self.gridSize-1)

       	self.domain  = np.array([
            [self.realRuntime[x] for x in 'xmin ymin zmin'.split()],
            [self.realRuntime[x] for x in 'xmax ymax zmax'.split()]
        ])

        self.domainSize = np.abs(self.domain[1]-self.domain[0])
        self.blockSize  = np.array([self.integerScalars[x] for x in 'nxb nyb nzb'.split()])
        self.cellSize   = self.domainSize / self.gridSize

    def get(self,dname):
        return self.h5file.get(dname).value
   
    def box(self, dname):
        # auxiliary variables for code clarity: shape: (3,)
        gridsize  = self.gridSize
        domsize   = self.domainSize
        offset    = -(self.domain[0])
        blksize   = self.blockSize

        # get block ids of desired refinement level
        rls  = self.get('refine level')
        rl   = rls.max()
        bids = [i for (i,x) in enumerate(rls) if x == rl] # block ids
        
        coords = (self.get('coordinates'))[bids] # shape: (#bids,3)
        blocks = (self.get(dname))[bids]         # shape: (#bids,nxb*nyb*nzb)

        # note: using numpy broadcasting kung-fu: shape: coords.shape
        positions = ((coords+offset)/domsize * gridsize).astype(np.int)

        box = np.zeros(gridsize)

        for bid, pos in enumerate(positions):
            I = np.array((pos-blksize//2,pos+blksize//2)).transpose()
            box[[slice(*i) for i in I]] = blocks[bid].transpose((2,1,0))

        return box 
    
    def list_datasets(self):
        return self.h5file.keys()

    def to_hdf(self,fpath):

        outfile = H5File(fpath,'w')

        outfile.create_dataset("DIMS", data=self.meta['grid size'])
        outfile.create_dataset("CELL VOL", data=self.meta['cell size'])
        outfile.create_dataset("BOX VOL", data=self.meta['domain size'])

        outfile.create_dataset("MIN BOUNDS", data=self.meta['min bounds'])
        outfile.create_dataset("MAX BOUNDS", data=self.meta['max bounds'])

        time = self.meta['real scalars']['time']
        step = float(self.meta['integer scalars']['nstep'])
        dt   = 0.0

        outfile.create_dataset("SIM_INFO", data=np.array([time,step,dt]))
        #GS = self.meta['grid size']
        #GV = GS[0]*GS[1]*GS[2]

        outfile.create_dataset("dens", data=self.get_box('dens').reshape(-1))
        outfile.create_dataset("velx", data=self.get_box('velx').reshape(-1))
        outfile.create_dataset("vely", data=self.get_box('vely').reshape(-1))
        outfile.create_dataset("velz", data=self.get_box('velz').reshape(-1))

        outfile.close()

# compatibility class for legacy code
class FlashFile(AMR):
    def get_box(self, *args, **kwargs):
        return self.box(*args, **kwargs) 

class HOPR:

    def get(self,dname):
        return self.h5file.get(dname).value
    
    def __init__(self,fpath):
        self.h5file = H5File(fpath,'r')
        self.elemInfo = self.get('ElemInfo')

        elemTypes = self.elemInfo[:,0].uniq
        if len(self.elemTypes) > 1:
            raise AssertionError('multiple element types detected')
        if self.elemTypes[0] != 108:
            raise AssertionError("type of all elements must be '108 aka. cube'")

        self.nodeCoords = self.get('NodeCoords')

       	self.domain  = np.array([
            [self.nodeCoords[i].min() for i in range(0,3)],
            [self.nodeCoords[i].max() for i in range(0,3)]
        ])

        self.domainSize = np.abs(self.domain[1]-self.domain[0])

    def cells(self):
        for element in enumerate(self.elemInfo):
            yield np.array(self.nodeCoords[el[4]], self.nodeCoords[el[5]-1])

class Flexi:

    def __init__(self, fpath, hopr):
        self.h5file = H5File(fpath,'r+')
        self.data = self.get('DG_Solution')
        self.hopr = hopr

    def save(self):
        pass

# hopr = Hopr(...)
# flexi = Flexi(hopr)
