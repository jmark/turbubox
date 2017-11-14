import numpy as np
import turbubox.h5 as h5
import turbubox.ulz as ulz
import turbubox.interpolate as interpolate
import turbubox.gausslobatto as gausslobatto

class BaseFile(h5.File):
    def __init__(self, fpath, mode='r', **kwargs):
        super().__init__(fpath, mode)

        self.framework  = 'couchdg'

        # support for older deprecated format
        if 'meta_num' in self.keys():
            self.meta = dict(
                tuple((k.strip().decode(),v.strip().decode()) for (k,v) in 
                    zip(self.get('meta_txt').attrs.get('keys'), self.get('meta_txt'))) + \
                tuple((k.strip().decode(),v) for (k,v) in 
                    zip(self.get('meta_num').attrs.get('keys'), self.get('meta_num')))
            )
        else:
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

class StructuredMeshFile(BaseFile):
    def __init__(self, fpath, mode='r'):
        super().__init__(fpath, mode)

        vs = np.array([self.get('/patches/'+pid+'/vertices')
                            for pid in self.file['patches'].keys()])
        self.domain = np.array(((np.min(vs[:,:,0]),np.min(vs[:,:,1])),
                                (np.max(vs[:,:,0]),np.max(vs[:,:,1]))))

        self.domsize = np.abs(self.domain[1]-self.domain[0])

    def as_box(self, ivar, Nvisu=None):
        return self.stitch(ivar, Nvisu)

    # depcrecated stitching routine
    def stitch_old(self, ivar, Nvisu=None, dname='state'):
        retv = None

        for pid in sorted(self.file['patches'].keys()):
            patch = self.get('/patches/'+pid+'/'+dname)

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
       
    def stitch_structured(self, ivar, Nvisu=None, dname='state'):
        NX_PATCHES = self.meta['mesh_nx_patches']
        NY_PATCHES = self.meta['mesh_ny_patches']

        Np = int(self.npoly)
        Nv = Nvisu if Nvisu else Np + 1

        xs = gausslobatto.mk_nodes(Np, self.nodetype)
        Xs = ulz.mk_body_centered_linspace(-1,1, Nv)

        Nx,Ny = self.meta['mesh_nx_cells'], self.meta['mesh_ny_cells'] 

        cloth = None
        for patchJ in range(NY_PATCHES):
            column = None
            for patchI in range(NX_PATCHES):
                patchid = patchJ * NX_PATCHES + patchI
                patch   = self.get('/patches/{:04d}/{}'.format(patchid,dname))
                patch   = patch[:,:,ivar,:,:].transpose(1,0,3,2).reshape((-1,Np+1,Np+1))
                patch   = interpolate.change_grid_space_2d(patch,xs,Xs).reshape(Nx,Ny,Nv,Nv)
                patch   = np.concatenate([np.concatenate(row,axis=1) for row in patch])
                column  = patch if column is None else np.concatenate((column,patch),axis=0)
            cloth = column if cloth is None else np.concatenate((cloth,column),axis=1)

        return cloth.T

    def stitch(self, *args,**kwargs):
        if 'mesh_type' in self.meta: return self.stitch_structured(*args,**kwargs)
        return self.stitch_old(*args,**kwargs)

    def get_prims(self, Nvisu=None, cons2prim=ulz.navier_conservative_to_primitive, gamma=5/3):
        cons = [None]*5
        cons[0] = self.stitch(0)
        cons[1] = self.stitch(1)
        cons[2] = self.stitch(2)
        cons[3] = 0.0
        cons[4] = self.stitch(3)

        return cons2prim(cons, gamma) 

Ribbon = StructuredMeshFile # support legacy code

# entry point and dispatch class
class File(BaseFile):
    def __init__(self, fpath, **kwargs):
        super().__init__(fpath, mode='r')
        if 'mesh_type' not in self.meta: # legacy support
            return StructuredMeshFile(fpath, **kwargs)
        elif self.meta['mesh_type'] == 'structured':
            return StructuredMeshFile(fpath, **kwargs)
        raise TypeError('Unsupported mesh type: {}'.format(self.meta['mesh_type']))
