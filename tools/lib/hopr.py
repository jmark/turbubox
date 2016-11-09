from h5 import H5File
import numpy as np

class MeshFile:

    def __init__(self,fpath, mode='r'):
        self.h5file     = H5File(fpath,mode)
        self.elemInfo   = self.get('ElemInfo')
        self.nodeCoords = self.get('NodeCoords')

       	self.domain  = np.array([
            [self.nodeCoords[:,i].min() for i in range(0,3)],
            [self.nodeCoords[:,i].max() for i in range(0,3)]
        ])

        self.domainsize = np.abs(self.domain[1]-self.domain[0])

    def get(self,dname):
        return self.h5file.get(dname)
    

class CartesianMeshFile(MeshFile):

    def __init__(self,fpath, mode='r'):
        super().__init__(fpath, mode)

        self.elemTypes = np.unique(self.elemInfo[:,0])
        if len(self.elemTypes) > 1:
            raise AssertionError('multiple element types detected: %s' % self.elemTypes)
        if self.elemTypes[0] != 108:
            raise AssertionError("type of all elements must be '108 aka. cube'")

        self.cellsize = np.abs(self.nodeCoords[7]-self.nodeCoords[0])
        self.gridsize = self.domainsize // self.cellsize
        self.nrelems  = len(self.elemInfo)
 
    def get_cell_coords(self):
        return (self.nodeCoords[:-7:8,:], self.nodeCoords[7::8,:])
