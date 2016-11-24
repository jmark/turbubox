from h5 import H5File
import numpy as np
import ulz
import interpolate
import gausslobatto

class File:
    def __init__(self, fpath, mesh, mode='r+'):
        self.h5file = H5File(fpath, mode)
        self.data = self.h5file.get('DG_Solution')
        self.mesh = mesh
        self.hopr = mesh
        self.attr = self.h5file.get("/").attrs
        self.nodetype   = self.attr['NodeType'][0].decode('utf-8').lower()
        self.Npoly      = self.attr['N'][0]
        self.Nout       = len(self.data[0,:,0,0,0])
        self.time       = self.attr['Time'][0]
        
        #self.params = dict((k,ulz.coerce(v)) for k,v in 
        #    [x.decode('utf8').split('=') for x in self.attr['Parameters']])

    # def flexi_to_box(self, iVar, Nvisu=None):
    #     if Nvisu is None:
    #         Nvisu = self.Nout

    #     xs = gausslobatto.mk_nodes(self.Nout-1,self.nodetype)
    #     Xs = ulz.mk_body_centered_linspace(-1,1,Nvisu)

    #     return interpolate.flexi_to_box(xs,Xs,self.data[:,:,:,:,iVar],self)

    def as_box(self, iVar, Nvisu=None):
        if Nvisu is None:
            Nvisu = self.Nout

        xs = gausslobatto.mk_nodes(self.Nout-1,self.nodetype)
        Xs = ulz.mk_body_centered_linspace(-1,1,Nvisu)

        elems = interpolate.change_grid_space(self.data[:,:,:,:,iVar].transpose(0,3,2,1),xs,Xs)
        return interpolate.elements_to_box(elems, self.mesh)

    def flexi_bo_box(self, iVar, Nvisu=None):
        return self.as_box(iVar, Nvisu)

    # def box_to_flexi(self, box, withBoundaryNodes=False):
    #     # init grid spaces
    #     xs  = ulz.mk_body_centered_linspace(-1,1,self.Nout, withBoundaryNodes)
    #     Xs  = gausslobatto.mk_nodes(self.Nout-1, self.nodetype)

    #     if withBoundaryNodes:
    #         box = ulz.wrap_in_guard_cells(box)        

    #     return interpolate.box_to_flexi(xs, Xs, box, self)

    def close(self):
        self.h5file.close()

    # provide context manager interface
    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()
        if isinstance(value,Exception):
            raise
