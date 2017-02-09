from h5 import H5File
import numpy as np
import ulz
import interpolate
import gausslobatto
import hopr
import os

def PeriodicBox(srcfp, meshfp=None, mode='r'):
    if meshfp is None:
        h5file = H5File(srcfp, 'r')
        meshfp = h5file.attrs['MeshFile'][0].decode('utf-8')
        h5file.close()
        oldcwd = os.path.realpath(os.curdir)
        #print(os.path.dirname(os.path.realpath(srcfp)))
        os.chdir(os.path.dirname(os.path.realpath(srcfp)))
        meshfp = os.path.realpath(meshfp)
        os.chdir(oldcwd)

    return File(srcfp, hopr.CartesianMeshFile(meshfp), mode)

class File:
    def __init__(self, fpath, mesh, mode='r'):
        self.h5file = H5File(fpath, mode)
        self.data = self.h5file.get('DG_Solution')
        self.mesh = mesh
        self.hopr = mesh
        self.attr = self.h5file.get("/").attrs
        self.nodetype   = self.attr['NodeType'][0].decode('utf-8').lower()
        self.Npoly      = self.attr['N'][0]
        self.Nout       = len(self.data[0,:,0,0,0])
        self.time       = self.attr['Time'][0]

        self.varnames   = 'dens momx momy momz eint'.split()

        self.domain     = self.mesh.domain
        self.domainsize = self.mesh.domainsize 
        self.cellsize   = self.mesh.cellsize
        self.cellvolume = np.prod(self.cellsize)
        
        #self.params = dict((k,ulz.coerce(v)) for k,v in 
        #    [x.decode('utf8').split('=') for x in self.attr['Parameters']])

    # def flexi_to_box(self, iVar, Nvisu=None):
    #     if Nvisu is None:
    #         Nvisu = self.Nout

    #     xs = gausslobatto.mk_nodes(self.Nout-1,self.nodetype)
    #     Xs = ulz.mk_body_centered_linspace(-1,1,Nvisu)

    #     return interpolate.flexi_to_box(xs,Xs,self.data[:,:,:,:,iVar],self)

    #def data(self, varname, Nvisu=None):
    #    iVar = self.varnames.index(varname)
    #    return self.as_box(iVar, Nvisu)

    def as_box(self, iVar, Nvisu=None):
        if Nvisu is None:
            Nvisu = self.Nout

        xs = gausslobatto.mk_nodes(self.Nout-1, self.nodetype)
        Xs = ulz.mk_body_centered_linspace(-1,1, Nvisu)

        elems = interpolate.change_grid_space(self.data[:,:,:,:,iVar].transpose(0,3,2,1),xs,Xs)
        return interpolate.elements_to_box(elems, self.mesh)

    def as_box_fv(self, iVar, Nvisu=None):
        if Nvisu is None:
            Nvisu = self.Nout

        xs = gausslobatto.mk_nodes(self.Nout-1, self.nodetype)
        Xs = ulz.mk_body_centered_linspace(-1,1, Nvisu)

        elems = self.data[:,:,:,:,iVar].transpose(0,3,2,1)
        return interpolate.elements_to_box(elems, self.mesh)

    def flexi_to_box(self, iVar, Nvisu=None):
        return self.as_box(iVar, Nvisu)

    def get_cons(self, Nvisu=None):
        return [self.flexi_to_box(i, Nvisu) for i in range(0,len(self.varnames))]

    def get_prims(self, Nvisu=None):
        cons = [self.as_box(i, Nvisu) for i in range(0,len(self.varnames))]
        return ulz.navier_conservative_to_primitive(cons)


    def get_cons_fv(self, Nvisu=None):
        return [self.as_box_fv(i, Nvisu) for i in range(0,len(self.varnames))]

    def get_prims_fv(self, Nvisu=None):
        cons = [self.as_box_fv(i, Nvisu) for i in range(0,len(self.varnames))]
        return ulz.navier_conservative_to_primitive(cons)

    def close(self):
        self.h5file.close()

    # provide context manager interface
    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()
        if isinstance(value,Exception):
            raise
