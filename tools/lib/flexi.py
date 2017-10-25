import h5
import numpy as np
import ulz
import interpolate
import gausslobatto
import hopr
import os

class File:
    def __init__(self, fpath, mesh, mode='r'):
        self.h5file     = h5.File(fpath, mode)
        self.data       = self.h5file.get('DG_Solution')
        self.mesh       = self.hopr = mesh
        self.attr       = self.h5file.get("/").attrs
        self.nodetype   = self.attr['NodeType'][0].decode('utf-8').lower()
        self.Npoly      = self.attr['N'][0]
        self.Nout       = len(self.data[0,:,0,0,0])
        self.time       = self.attr['Time'][0]

        self.var2idx    = dict((k,v) for v,k in enumerate('dens momx momy momz ener magx magy magz'.split()))

        self.domain     = self.mesh.domain
        self.domainsize = self.mesh.domainsize 
        self.cellsize   = self.mesh.cellsize
        self.cellvolume = np.prod(self.cellsize)
        
    def as_box(self, iVar, Nvisu=None):
        if isinstance(iVar, str): iVar = self.var2idx[iVar]
        if Nvisu is None: Nvisu = self.Nout

        xs = gausslobatto.mk_nodes(self.Nout-1, self.nodetype)
        Xs = ulz.mk_body_centered_linspace(-1,1, Nvisu)

        elems = interpolate.change_grid_space(self.data[:,:,:,:,iVar].transpose(0,3,2,1),xs,Xs)
        return interpolate.elements_to_box(elems, self.mesh)

    def as_box_fv(self, iVar, Nvisu=None):
        """This routine works the same as 'as_box' but recognizes finite-volume cells and
           treats them appropiately."""

        if Nvisu is None:
            Nvisu = self.Nout

        # Which elements are FV cells?
        elemData = self.h5file.get('ElemData')
        FVs = elemData[:,2].astype(np.int32)
    
        xs = gausslobatto.mk_nodes(self.Nout-1, self.nodetype)
        Xs = ulz.mk_body_centered_linspace(-1,1, Nvisu)

        elems = interpolate.change_grid_space_dg_fv(self.data[:,:,:,:,iVar].transpose(0,3,2,1),xs,Xs,FVs)
        return interpolate.elements_to_box(elems, self.mesh)

    def get_prims(self, Nvisu=None, cons2prim=ulz.navier_conservative_to_primitive, gamma=5/3):
        if Nvisu is None:
            Nvisu = self.Nout

        # Which elements are FV cells?
        elemData = self.h5file.get('ElemData')
        FVs = elemData[:,2]
    
        xs = gausslobatto.mk_nodes(self.Nout-1, self.nodetype)
        Xs = ulz.mk_body_centered_linspace(-1,1, Nvisu)

        cons  = [self.data[:,:,:,:,i] for i in range(0,len(self.varnames))] 
        prims = cons2prim(cons, gamma)

        retval = []
        for iVar in range(0,len(self.varnames)):
            elems = interpolate.change_grid_space_dg_fv(prims[iVar].transpose(0,3,2,1),xs,Xs,FVs)
            retval.append(interpolate.elements_to_box(elems, self.mesh))

        return retval

    def get_cons(self, Nvisu=None):
        if Nvisu is None:
            Nvisu = self.Nout

        # Which elements are FV cells?
        elemData = self.h5file.get('ElemData')
        FVs = elemData[:,2]
    
        xs = ulz.mk_body_centered_linspace(-1,1, Nvisu)
        Xs = gausslobatto.mk_nodes(self.Nout-1, self.nodetype)

        cons  = [self.data[:,:,:,:,i] for i in range(0,len(self.varnames))] 

        retval = []
        for iVar in range(0,len(self.varnames)):
            elems = interpolate.change_grid_space_fv_dg(cons[iVar].transpose(0,3,2,1),xs,Xs,FVs)
            retval.append(interpolate.elements_to_box(elems, self.mesh))

        return retval

    def close(self):
        self.h5file.close()

    # provide context manager interface
    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()
        if isinstance(value,Exception):
            raise

    # convenience methods
    def flexi_to_box(self, iVar, Nvisu=None):
        return self.as_box(iVar, Nvisu)

    def get_cons(self, Nvisu=None):
        return [self.as_box_fv(i, Nvisu) for i in range(0,len(self.varnames))]
