from h5 import H5File
import numpy as np
import sys
import ulz

class File:
    def __init__(self,fpath, mode='r'):
        self.h5file = H5File(fpath, mode)

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
                for N,n in zip('nblockx nblocky nblockz'.split(), 'nxb nyb nzb'.split())]).astype(np.int)

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

    def close(self):
        self.h5file.close()

    def __delete__(self):
        self.close()

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

class FakeFile:
    """This class is a mess, but it works. No docs here. RTFSC."""
    def __init__(self, domain, boxdims, fillby='constant', blockdims=None):
       	self.domain     = np.array(domain)
        self.gridsize   = np.array(boxdims).astype(np.int)
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

        def exp3D(A,sigma,p,x,y,z):
            return A * np.exp(-((x-p[0])**2 + (y-p[1])**2 + (z-p[2])**2)/sigma/2.)

        X,Y,Z = np.meshgrid(*tuple(ulz.mk_body_centered_linspace(*d, s) for d,s in zip(dom, grd)), indexing='ij')

        if dname == 'dens':

            D1 =  1.0
            D2 =  2.0

            ret = 0.1 + 0.9 * np.exp(-((Y-0.5)**2)/2/0.01)
            #ret = D1 * np.ones_like(X)
            #ret[np.where((6/16 <= Y) * (Y <= 10/16))] = D2 
            #ret[np.where((6/16 <= Y) * (Y <= 10/16))] = D2 

            # --------------------------------------------------------------- #
            # if fillby == 'constant':
            #     ret[:] = 1.0

            # elif fillby == 'planeX':
            #     ret[:] = X

            # elif fillby == 'planeX+':
            #     ret = np.where(
            #         (2/8-1/grd[0]<X) * (X<6/8+1/grd[0]) * \
            #         (2/8-1/grd[1]<Y) * (Y<6/8+1/grd[1]) * \
            #         (2/8-1/grd[2]<Z) * (Z<6/8+1/grd[2]), X,0*X)

            # elif fillby == 'planeXYZ':
            #     ret = np.where(
            #         (1/8-1/grd[0]<X) * (X<7/8+1/grd[0]) * \
            #         (1/8-1/grd[1]<Y) * (Y<7/8+1/grd[1]) * \
            #         (1/8-1/grd[2]<Z) * (Z<7/8+1/grd[2]), X+Y+Z,0*X)

            # elif fillby == 'plane+wiggle':
            #     ret = np.where(
            #         (1/16-1/grd[0]<X) * (X<15/16+1/grd[0]) * \
            #         (1/16-1/grd[1]<Y) * (Y<15/16+1/grd[1]) * \
            #         (1/16-1/grd[2]<Z) * (Z<15/16+1/grd[2]),
            #         X+Y+Z + 0.5*self.wiggle(X+1/16,Y+1/16,Z+1/16),0*X)

            # elif fillby == 'gaussianXYZ':
            #     ret = np.exp(-((X-0.5)**2 + (Y-0.5)**2 + (Z-0.5)**2)/2/0.02)

            # elif fillby == 'stepsXYZ':
            #     ret = X + 100*X*(np.sin(4*2*np.pi*X) + 0.2*np.cos(20*2*np.pi*X) + 0.1*np.sin(20*2*np.pi*X))**2 

            # else:
            #     raise NotImplementedError('unknow fillby: %s' % fillby)

        elif dname == 'pres':
            # constant
            #ret[:] = 1.0
            # plane
            #ret = X+Y

            # --------------------------------------------------------------- #
            # shear flow
            ret = np.ones_like(X)
            #ret[np.where((0.375 <= Y) * (Y <= 0.625) * (0.375 <= Z) * (Z <= 0.625))] = 1.0
            #ret[np.where((6/16 <= Y) * (Y <= 10/16))] = 5.0
            #ret[np.where((Y < 5/16) + (11/16 < Y))] = -5.0

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

        elif dname == 'velx':
            ## smashing balls
            #V0 = -0.8
            #V1 = -2 * V0 
            #s1 = 0.01
            #p1 = [0.5,0.5,0.5] 

            #ret = V0 * np.ones_like(X) + exp3D(V1, s1, p1, X,Y,Z)

            # --------------------------------------------------------------- #
            ## beam
            #ret = 0.1 * (2 * np.random.rand(*X.shape) - 1)
            #ret[np.where((0.375 <= Y) * (Y <= 0.625) * (0.375 <= Z) * (Z <= 0.625))] = 1.0
            #ret[np.where((Y < 23/64) + (41/64 < Y))] = -5.0

            U1 = -0.2
            U2 =  0.4

            #ret = U1 * np.ones_like(X)
            #ret[np.where((6/16 <= Y) * (Y <= 10/16))] = U2 

            ret = U1 * np.ones_like(X) + U2 * np.exp(-((Y-0.5)**2)/2/0.01)

        elif dname == 'vely':
            U1 =  0.01
            ret = U1 * np.sin(4*np.pi*X)

        elif dname in 'velz magx magy magz'.split():
            ret[:] = 0.0

        else:
            raise KeyError('%s not found!' % dname) 

        return ret
