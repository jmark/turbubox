from h5 import H5File
import numpy as np
import sys
import ulz

class File:
    def __init__(self,fpath):
        self.h5file = H5File(fpath,'r')

        # transform meta information to python data types
        self.siminfo        = self.get('sim info')
        
        self.refine_levels  = self.get('refine level')
        self.maxrefinelevel = self.refine_levels.max()

        self.integerscalars = H5File.dataset_to_dict(self.get('integer scalars'))
        self.realscalars    = H5File.dataset_to_dict(self.get('real scalars'))
        self.integerruntime = H5File.dataset_to_dict(self.get('integer runtime parameters'))
        self.realruntime    = H5File.dataset_to_dict(self.get('real runtime parameters'))

        # figure out global grid size
        self.gridsize = np.array([
            self.integerruntime[N] * self.integerscalars[n]*2**(self.maxrefinelevel-1) 
                for N,n in zip('nblockx nblocky nblockz'.split(), 'nxb nyb nzb'.split())])

        # handle uniform grid case
        if str(self.siminfo['setup call'][0]).find('+ug') >= 0:
            self.gridsize *= np.array([self.integerscalars[key] for key in 'iprocs jprocs kprocs'.split()])

        self.grid = np.array([[0,0,0], self.gridsize-1])

       	self.domain  = np.array([
            [self.realruntime[x] for x in 'xmin ymin zmin'.split()],
            [self.realruntime[x] for x in 'xmax ymax zmax'.split()]
        ])

        self.domainsize = np.abs(self.domain[1]-self.domain[0])
        self.blocksize  = np.array([self.integerscalars[x] for x in 'nxb nyb nzb'.split()])
        self.cellsize   = self.domainsize / self.gridsize

    def get(self,dname):
        return self.h5file.get(dname).value
   
    def data(self, dname):
        return self.get_data(dname)

    def get_data(self, dname):
        # auxiliary variables for code clarity: shape: (3,)
        gridsize  = self.gridsize
        domsize   = self.domainsize
        offset    = -(self.domain[0])
        blksize   = self.blocksize

        # get block ids of desired refinement level
        rls  = self.get('refine level')
        rl   = rls.max()
        bids = [i for (i,x) in enumerate(rls) if x == rl] # filter desired blocks
        
        coords = (self.get('coordinates'))[bids] # shape: (#bids,3)
        blocks = (self.get(dname))[bids]         # shape: (#bids,nxb*nyb*nzb)

        # note: using numpy broadcasting kung-fu: shape: coords.shape
        positions = np.round((coords+offset)/domsize * gridsize).astype(np.int)

        box = np.zeros(gridsize)
        for bid, pos in enumerate(positions):
            I = np.array((pos-blksize//2,pos+blksize//2)).transpose()
            box[[slice(*i) for i in I]] = blocks[bid].transpose((2,1,0))

        return box 

    def set_data(self, dname, fun):
        # auxiliary variables for code clarity: shape: (3,)
        gridsize  = self.gridsize
        domsize   = self.domainsize
        offset    = -(self.domain[0])
        blksize   = self.blocksize

        linspace  = ulz.mk_body_centered_linspace
        X,Y,Z     = np.meshgrid(*tuple(linspace(-1, 1, nb) for nb in blksize))

        coords = self.get('coordinates')    # shape: (#bids,3)
        blocks = self.get(dname)            # shape: (#bids,nxb*nyb*nzb)

        for bid, coord in enumerate(coords):
            mg = np.meshgrid(*tuple(linspace(cs[0], cs[1], nb) for cs,nb in zip(coord,blksize)), indexing='ij')
            blocks[bid] = fun(*mg).transpose((2,1,0))

    def list_datasets(self):
        return self.h5file.keys()

    def to_hdf(self,fpath): # WIP!

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
        outfile.create_dataset("dens", data=self.get_box('dens').reshape(-1))
        outfile.create_dataset("velx", data=self.get_box('velx').reshape(-1))
        outfile.create_dataset("vely", data=self.get_box('vely').reshape(-1))
        outfile.create_dataset("velz", data=self.get_box('velz').reshape(-1))

        outfile.close()
