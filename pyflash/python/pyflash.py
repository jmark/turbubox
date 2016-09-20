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
            
class FlashFile:

    def __init__(self,fpath):
        self.h5file = H5File(fpath,'r')
        self.meta = {}
        self.meta['max refine level'] = self.get_dataset('refine level').max()
        self.meta['integer scalars']  = self.list_to_dict('integer scalars')
        self.meta['real scalars']     = self.list_to_dict('real scalars')
        self.meta['integer runtime']  = self.list_to_dict('integer runtime parameters')
        self.meta['real runtime']     = self.list_to_dict('real runtime parameters')
        self.meta['sim info']         = self.get_dataset('sim info')
   
        self.meta['grid size'] = np.array([
            self.meta['integer scalars'][x]*2**(self.meta['max refine level']-1)
                for x in 'nxb nyb nzb'.split()])

        if str(self.meta['sim info']['setup call'][0]).find('+ug') >= 0:
            self.meta['grid size'] *= 2 

       	self.meta['min bounds']  = np.array([self.meta['real runtime'][x] for x in 'xmin ymin zmin'.split()])
        self.meta['max bounds']  = np.array([self.meta['real runtime'][x] for x in 'xmax ymax zmax'.split()])
        self.meta['domain size'] = np.abs(self.meta['min bounds']-self.meta['max bounds'])
        self.meta['block size']  = np.array([self.meta['integer scalars'][x] for x in 'nxb nyb nzb'.split()])
        self.meta['cell size']   = self.meta['domain size'] / self.meta['grid size']

    def get_dataset(self,dname):
        return self.h5file.get(dname).value
    
    def list_to_dict(self,dname):
        return dict([(k.strip().decode(),v) for (k,v) in self.h5file.get(dname)])
    
    def get_box(self,dname):
        # auxiliary variables for code clarity: shape: 3,
        gridsize  = self.meta['grid size']
        domsize   = self.meta['domain size']
        offset    = -(self.meta['min bounds'])
        blksize   = self.meta['block size']

        # get block ids of desired refinement level
        rl   = self.meta['max refine level'] 
        rls  = self.get_dataset('refine level')
        bids = [i for (i,x) in enumerate(rls) if x == rl]
        
        coords = (self.get_dataset('coordinates'))[bids] # shape: #bids,3
        blocks = (self.get_dataset(dname))[bids] # shape: #bids,nxb*nyb*nzb

        # note: using kung-fu broadcasting rules of numpy
        positions = ((coords + offset)/domsize * gridsize).astype(np.int)

        BOX = np.zeros(gridsize)

        for bid, pos in enumerate(positions):
            (i_min,i_max) = (pos[0] - blksize[0]//2, pos[0] + blksize[0]//2) 
            (j_min,j_max) = (pos[1] - blksize[1]//2, pos[1] + blksize[1]//2) 
            (k_min,k_max) = (pos[2] - blksize[2]//2, pos[2] + blksize[2]//2) 

            #print(i_min,i_max,j_min,j_max,k_min,k_max)

            BOX[i_min:i_max,j_min:j_max,k_min:k_max] = blocks[bid].transpose((2,1,0))

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
