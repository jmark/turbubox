#!/usr/bin/env pyturbubox

import numpy as np
import h5 as h5
import ulz as ulz
import interpolate as interpolate
import gausslobatto as gausslobatto

class BaseFile(h5.File):
    def __init__(self, fpath, mode='r', **kwargs):
        super().__init__(fpath, mode)

        self.framework  = 'couchdg'

        # support older, soon to be deprecated, format
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

        self.profiles = dict(
            tuple((k.strip().decode(),v-1) for (k,v) in 
                zip(self.get('profile_indices').attrs.get('keys'), self.get('profile_indices')))
        )

    def topology(self):
        for pid in sorted(self.file['patches'].keys()):
            yield self.get('/patches/'+pid+'/topology')

    def vertices(self):
        for pid in sorted(self.file['patches'].keys()):
            yield self.get('/patches/'+pid+'/vertices')

    def patches(self, ivar, Nprofilevisu=None):
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
       
    def stitch_structured_old(self, ivar, Nvisu=None, dname='state', imspace=None):
        NX_PATCHES = self.meta['mesh_nx_patches']
        NY_PATCHES = self.meta['mesh_ny_patches']

        Np = int(self.npoly)
        Nv = Nvisu if Nvisu else Np + 1

        xs = gausslobatto.mk_nodes(Np, self.nodetype)
        #xs = np.array([-0.875, -0.5  , -0.125,  0.125,  0.5  ,  0.875])
        if imspace is None:
            Xs = ulz.mk_body_centered_linspace(-1,1, Nv)
        else:
            Xs = imspace

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

    def stitch_structured(self, ivar, Nvisu=None, mpoly=None, dname='state', imspace=None):
        NX_PATCHES = self.meta['mesh_nx_patches']
        NY_PATCHES = self.meta['mesh_ny_patches']

        npoly = self.npoly
        mpoly = mpoly if mpoly else self.mpoly

        Nv = Nvisu if Nvisu else 2*(npoly + 1)
        Nw = Nv

        nn = int((npoly+1)/(mpoly+1))
        ns = np.array(gausslobatto.mk_nodes(mpoly,ntype=self.nodetype))
        xs = ulz.mk_body_centered_linspace(-1,1,Nv)
        ys = np.linspace(-1,1,nn,endpoint=False)
        Xs = np.array([-1 + nn*abs(ys[int(np.floor(0.5*nn*(y+1)))] - y) for y in xs])
        ks = list(range(mpoly+1)) * nn

        def krondeltas(I,i):
            j = int(np.floor(0.5*nn*(xs[I]+1)))
            l = list(range((mpoly+1)*j,(mpoly+1)*j+mpoly+1))
            return int(i in l)

        LL = np.array([[[[
                krondeltas(I,i)*krondeltas(J,j)*\
                gausslobatto.LagrangePolynomial(ns,ks[i],Xs[I])*\
                gausslobatto.LagrangePolynomial(ns,ks[j],Xs[J])
                    for j in range(npoly+1)]
                    for i in range(npoly+1)]
                    for J in range(len(Xs))]
                    for I in range(len(Xs))])

        Np    = int(self.npoly)
        Nx,Ny = self.meta['mesh_nx_cells'], self.meta['mesh_ny_cells'] 
        Nv,Nw = LL.shape[0],LL.shape[1]

        cloth = None
        for patchJ in range(NY_PATCHES):
            column = None
            for patchI in range(NX_PATCHES):
                patchid = patchJ * NX_PATCHES + patchI
                patch   = self.get('/patches/{:04d}/{}'.format(patchid,dname))
                patch   = patch[:,:,ivar,:,:].transpose(1,0,3,2).reshape((-1,Np+1,Np+1))
                patch   = interpolate.change_grid_space_2d_2(LL,patch)
                patch   = patch.reshape(Nx,Ny,Nv,Nw)
                patch   = np.concatenate([np.concatenate(row,axis=1) for row in patch])
                column  = patch if column is None else np.concatenate((column,patch),axis=0)
            cloth = column if cloth is None else np.concatenate((cloth,column),axis=1)

        return cloth.T

    def stitch(self, *args,**kwargs):
        if 'mesh_type' in self.meta: return self.stitch_structured(*args,**kwargs)
        return self.stitch_old(*args,**kwargs)

    def get_prims(self, Nvisu=None, imspace=None, cons2prim=ulz.navier_conservative_to_primitive, gamma=5/3):
        cons = [None]*5
        cons[0] = self.stitch(0, Nvisu=Nvisu, imspace=imspace)
        cons[1] = self.stitch(1, Nvisu=Nvisu, imspace=imspace)
        cons[2] = self.stitch(2, Nvisu=Nvisu, imspace=imspace)
        cons[3] = 0.0
        cons[4] = self.stitch(3, Nvisu=Nvisu, imspace=imspace)

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

if __name__ == '__main__':
    import sys
    fp = sys.argv[1]
    fh = StructuredMeshFile(fp)

    mpoly = 3
    Nvisu = 17

    data = fh.stitch_structured_2(0,Nvisu=Nvisu,mpoly=mpoly,dname='state')
    #data = fh.stitch_structured(0,Nvisu=24,dname='state')

    import matplotlib
    matplotlib.rcParams.update({'font.size': 9})
    import matplotlib.pyplot as plt

    #fig = plt.figure(1,figsize=(geo[0]/dpi, geo[1]/dpi), dpi=dpi)
    fig  = plt.figure()

    plt.title('density: t = {:.2f}'.format(fh.time))
    plt.imshow(
        data,
        #extent = extent,
        origin='lower',
        interpolation = None,
        cmap = plt.get_cmap('cubehelix'),
    )
    plt.colorbar()

    plt.show()
