import sys
import ulz
import numpy as np

import h5 as h5
import interpolate as itpl

class File(h5.File):
    def __init__(self, fpath, mode='r', **kwargs):
        super().__init__(fpath, mode)
        self.framework      = 'flash'

        # transform meta information to python data types
        self.siminfo        = self.get('sim info')
        
        self.refine_levels  = self.get('refine level')[()]
        self.maxrefinelevel = np.max(self.refine_levels)
        self.is_multilevel  = any(level > 1 for level in self.refine_levels)

        self.realscalars    = h5.dataset_to_dict(self.get('real scalars'))
        self.realruntime    = h5.dataset_to_dict(self.get('real runtime parameters'))
        self.integerscalars = h5.dataset_to_dict(self.get('integer scalars'))
        self.integerruntime = h5.dataset_to_dict(self.get('integer runtime parameters'))

        for key in 'nblockx nblocky nblockz'.split():
            if not key in self.integerruntime:
                self.integerruntime[key] = 1

        if self.integerscalars['nzb'] == 1:
            self.ndims = 2
        else:
            self.ndims = 3

        self.coords         = self.get('coordinates')[()]
        self.gridsize       = self.calc_gridsize(self.maxrefinelevel)
        self.grid           = np.array([[0,0,0], self.gridsize-1])

        self.domain  = np.array([
            [self.realruntime[x] for x in 'xmin ymin zmin'.split()],
            [self.realruntime[x] for x in 'xmax ymax zmax'.split()]
        ])

        self.domainsize = np.abs(self.domain[1]-self.domain[0])
        self.domsize    = self.domainsize
        self.blocksize  = np.array([self.integerscalars[x] for x in 'nxb nyb nzb'.split()])
        self.cellsize   = self.domainsize / self.gridsize
        self.cellvolume = np.prod(self.cellsize)

        # shortcut to general parameters useful in analysis
        self.params = dict()
        self.params['time']  = self.realscalars['time']
        self.params['dt']    = self.realscalars['dt']
        self.params['gamma'] = self.realruntime['gamma']
        self.params['kappa'] = self.realruntime['gamma'] # synonym
        self.gamma           = self.params['gamma']
        self.time            = self.params['time']

        if self.ndims == 2:
            self.extent = tuple(self.realruntime[k] for k in 'xmin xmax ymin ymax'.split())
        else:
            self.extent = tuple(self.realruntime[k] for k in 'xmin xmax ymin ymax zmin zmax'.split())

    def __getattr__(self, name):
        if name in self.realscalars:
            return self.realscalars[name]
        if name in self.realruntime:
            return self.realruntime[name]
        if name in self.integerscalars:
            return self.integerscalars[name]
        if name in self.integerruntime:
            return self.integerruntime[name]

        raise AttributeError('Unknown attritube: {}'.format(name))

    def data(self,dname):
        return self.get_data(dname)

    def calc_gridsize(self,rlevel):
        gridsize = np.array([self.integerruntime[N] * self.integerscalars[n]*2**(rlevel-1) 
                for N,n in zip('nblockx nblocky nblockz'.split(), 'nxb nyb nzb'.split())]).astype(np.int)

        if self.ndims == 2:
            gridsize = np.array((gridsize[0],gridsize[1],1))

        if not self.is_multilevel: # handle uniform grid
            gridsize *= np.array([self.integerscalars[key] for key in 'iprocs jprocs kprocs'.split()])

        return gridsize

    def get_data(self,dname,shape=None,method='nearest'):
        levels = self.get('refine level')[()]
        coords = self.get('coordinates')[()]
        bsizes = self.get('block size')[()]
        ntype  = self.get('node type')[()]

        domsize = self.domainsize
        if self.ndims == 2:
            domsize[2] = 1.0

        coords -= self.domain[0]
        coords /= domsize
        bsizes /= domsize

        if shape == None:
            shape = 2**(self.maxrefinelevel-1) * self.blocksize

        if self.ndims == 2:
            image = np.zeros(shape[0:2])
            blocks = np.transpose(self.get(dname),(0,3,2,1))
            itpl.cells_to_image_2d(ntype,coords,bsizes,blocks,image,method=method)

        return image

    def as_box(self,dname):
        return self.get_data(dname)

    def get_prims(self,**kwargs):
        return [self.get_data(dname, **kwargs) for dname in 'dens velx vely velz pres'.split()]

    def get_cons(self,gamma=5./3.):
         return ulz.navier_primitive_to_conservative(self.get_prims(), gamma)

    # experimental!
    def set_data(self,dname,box):
        if self.is_multilevel:
            raise NotImplementedError('Setting data for multilevel grids (AMR) is not supported yet!')

        if box.shape != tuple(self.gridsize):
            raise ValueError('Given box shape does not match gridsize!')

        # auxiliary variables for code clarity: shape: (3,)
        gridsize  = self.gridsize
        domsize   = self.domainsize
        offset    = -(self.domain[0])
        blksize   = self.blocksize

        linspace  = ulz.mk_body_centered_linspace
        X,Y,Z     = np.meshgrid(*tuple(linspace(-1, 1, nb) for nb in blksize))

        coords = self.get('coordinates')    # shape: (#bids,3)
        blocks = self.h5file.get(dname)     # shape: (#bids,nxb*nyb*nzb)

        # note: using numpy broadcasting kung-fu: shape: coords.shape
        positions = np.round((coords+offset)/domsize * gridsize).astype(np.int)

        for bid, pos in enumerate(positions):
            I = np.array((pos-blksize//2,pos+blksize//2)).transpose()
            blocks[bid] = box[[slice(*i) for i in I]].transpose((2,1,0))

        #blocks.file.flush()

    def list_datasets(self):
        return self.h5file.keys()

    def to_hdf(self,fpath): # WIP!
        outfile = h5.File(fpath,'w')

        outfile.create_dataset("DIMS", data=self.meta['grid size'])
        outfile.create_dataset("CELL VOL", data=self.meta['cell size'])
        outfile.create_dataset("BOX VOL", data=self.meta['domain size'])

        outfile.create_dataset("MIN BOUNDS", data=self.meta['min bounds'])
        outfile.create_dataset("MAX BOUNDS", data=self.meta['max bounds'])

        time = self.meta['real scalars']['time']
        step = float(self.meta['integer scalars']['nstep'])
        dt   = 0.0

        outfile.create_dataset("SIM_INFO", data=np.array([time,step,dt]))
        outfile.create_dataset("dens", data=self.get_box('dens').reshape(-1))
        outfile.create_dataset("velx", data=self.get_box('velx').reshape(-1))
        outfile.create_dataset("vely", data=self.get_box('vely').reshape(-1))
        outfile.create_dataset("velz", data=self.get_box('velz').reshape(-1))

        outfile.close()

if __name__ == '__main__':
    fp = sys.argv[1]
    fh = File(fp)

    # dens = fh.get_data('dens')
    # print(dens.shape)

    levels = fh.get('refine level')[()]
    coords = fh.get('coordinates')[()]
    bsizes = fh.get('block size')[()]
    ntype  = fh.get('node type')[()]

    domsize = fh.domainsize
    coords -= fh.domain[0]
    coords /= domsize
    bsizes /= domsize

    p = np.array([0.0,0.0,0.1])
    u = np.array([1.0,0.0,0.0])
    v = np.array([0.0,1.0,0.0])

    nedges = itpl.plane_morton_to_coords(ntype,coords,bsizes, p,u,v, edges=None)
    edges = np.zeros([nedges,2,3])
    nedges = itpl.plane_morton_to_coords(ntype,coords,bsizes, p,u,v, edges=edges)

    print(nedges)
    print(edges)
