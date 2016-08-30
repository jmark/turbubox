import h5py as h5
import numpy as np
import sys

class H5File(h5.File):

    # def __init__(self,fpath,mode):
    #     try:
    #         super().__init__(fpath,mode)
    #     except OSError as e:
    #         print(sys.stderr,"OSError: ",e)
    #         sys.exit(1)
    
    def get(self,dname):
        db = super(H5File,self).get(dname)
        if db:
            return db
        else:
            raise KeyError(dname)
            
def distance(a,b):
    return np.abs(a-b)
            
class FlashFile:

    def __init__(self,fpath):
        self.h5file = H5File(fpath,'r')
        self.meta = {}
        self.meta['max refine level'] = self.get_dataset('refine level').max()
        self.meta['integer scalars']  = self.list_to_dict('integer scalars')
        self.meta['real scalars']     = self.list_to_dict('real scalars')
        self.meta['integer runtime']  = self.list_to_dict('integer runtime parameters')
        self.meta['real runtime']     = self.list_to_dict('real runtime parameters')
        
        self.meta['grid size'] = np.array([
            self.meta['integer scalars']['nxb']*2**(self.meta['max refine level']-1),
            self.meta['integer scalars']['nyb']*2**(self.meta['max refine level']-1),
            self.meta['integer scalars']['nzb']*2**(self.meta['max refine level']-1)
        ])
        
       	self.meta['min bounds'] = np.array([
            self.meta['real runtime']['xmin'],
            self.meta['real runtime']['ymin'],
            self.meta['real runtime']['zmin']
        ])

        self.meta['max bounds'] = np.array([
            self.meta['real runtime']['xmax'],
            self.meta['real runtime']['ymax'],
            self.meta['real runtime']['zmax']
        ])

        self.meta['domain size'] = distance(self.meta['min bounds'],self.meta['max bounds'])
        self.meta['cell size']   = self.meta['domain size'] / self.meta['grid size']

    def get_dataset(self,dname):
        return self.h5file.get(dname).value
    
    def list_to_dict(self,dname):
        return dict([(k.strip().decode(),v) for (k,v) in self.h5file.get(dname)])
    
    def refine_level_indices(self):
        return [i for (i,v) in enumerate(self.get_dataset('refine level')) if v == self.meta['max refine level']]
    
    def block_index_to_range(self,bix,dim):
        n_b = [self.meta['integer scalars'][el] for el in ['nxb', 'nyb', 'nzb']]
        #return (bix[dim] * n_b[dim], (bix[dim] + 1)*n_b[dim])
        return (bix[dim], bix[dim] + n_b[dim])
        
    def get_box(self,dname):
        BOX = np.zeros(tuple(self.meta['grid size'][i] for i in range(0,3)))
        
        rli = self.refine_level_indices()
        
        # following is a mess like hell ...
        coords    = (self.get_dataset('coordinates'))[rli]
        #subdomsize   = (self.get_dataset('block size'))[rli]
        boxsize   = self.meta['grid size'][0]
        blksize   = self.meta['integer scalars']['nxb']
        domsize   = self.meta['domain size'][0]
        offset    = (0.0 - self.meta['real runtime']['xmin']) * boxsize/domsize
        bindices  = ((coords * boxsize/domsize - blksize/2 + offset).astype(np.int))
        blocks    = (self.get_dataset(dname))[rli]           
        
        for bix,bindex in enumerate(bindices):
            i_min,i_max = self.block_index_to_range(bindex,0)
            j_min,j_max = self.block_index_to_range(bindex,1)
            k_min,k_max = self.block_index_to_range(bindex,2)
            
            BOX[i_min:i_max,j_min:j_max,k_min:k_max] = blocks[bix].transpose((2,1,0))
            
        return BOX
    
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
