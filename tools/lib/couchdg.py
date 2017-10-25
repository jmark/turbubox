import h5
import ulz
import interpolate
import gausslobatto

import numpy as np

class File(h5.File):
    def __init__(self, fpath, mode='r'):
        super().__init__(fpath, mode)

        self.meta = dict(
            tuple((k.strip().decode(),v.strip().decode()) for (k,v) in 
                zip(self.get('meta_txt').attrs.get('keys'), self.get('meta_txt'))) + \
            tuple((k.strip().decode(),v) for (k,v) in 
                zip(self.get('meta_int').attrs.get('keys'), self.get('meta_int'))) + \
            tuple((k.strip().decode(),v) for (k,v) in 
                zip(self.get('meta_flt').attrs.get('keys'), self.get('meta_flt')))
        )

        for k,v in self.meta.items(): setattr(self, k, v)

    def topology(self):
        for pid in sorted(self.file['patches'].keys()):
            yield self.get('/patches/'+pid+'/topology')

    def vertices(self):
        for pid in sorted(self.file['patches'].keys()):
            yield self.get('/patches/'+pid+'/vertices')

    def patches(self, ivar, Nvisu=None):
        for pid in sorted(self.file['patches'].keys()):
            patch = self.get('/patches/'+pid+'/state')

            Np = int(self.npoly)
            Nv = Nvisu if Nvisu else Np + 1
            Ny,Nx = patch.shape[0:2]

            xs = gausslobatto.mk_nodes(Np, self.nodetype)
            Xs = ulz.mk_body_centered_linspace(-1,1, Nv)

            temp = patch[:,:,ivar,:,:].transpose(1,0,3,2).reshape((-1,Np+1,Np+1))
            temp = interpolate.change_grid_space_2d(temp,xs,Xs).reshape(Nx,Ny,Nv,Nv)
            yield np.concatenate([np.concatenate(row,axis=1) for row in temp]).T

class Ribbon(File):
    def __init__(self, fpath, mode='r'):
        super().__init__(fpath, mode)

        vs = np.array([self.get('/patches/'+pid+'/vertices')
                            for pid in self.file['patches'].keys()])
        self.domain = np.array(((np.min(vs[:,:,0]),np.min(vs[:,:,1])),
                                (np.max(vs[:,:,0]),np.max(vs[:,:,1]))))

        self.domsize = np.abs(self.domain[1]-self.domain[0])

    def stitch(self, ivar, Nvisu=None):
        retv = None

        for pid in sorted(self.file['patches'].keys()):
            patch = self.get('/patches/'+pid+'/state')

            Np = int(self.npoly)
            Nv = Nvisu if Nvisu else Np + 1
            Ny,Nx = patch.shape[0:2]

            xs = gausslobatto.mk_nodes(Np, self.nodetype)
            Xs = ulz.mk_body_centered_linspace(-1,1, Nv)

            temp = patch[:,:,ivar,:,:].transpose(1,0,3,2).reshape((-1,Np+1,Np+1))
            temp = interpolate.change_grid_space_2d(temp,xs,Xs).reshape(Nx,Ny,Nv,Nv)
            temp = np.concatenate([np.concatenate(row,axis=1) for row in temp])

            retv = temp if retv is None else np.concatenate((retv,temp),axis=0)

        return retv.T
