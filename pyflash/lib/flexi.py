from h5 import H5File
import numpy as np

class File:
    def __init__(self, fpath, mesh, mode='r+'):
        self.h5file = H5File(fpath, mode)
        self.data = self.get('DG_Solution')
        self.mesh = mesh
        self.attr = self.h5file.get("/").attrs
        self.N = self.attr['N'][0]

    def get(self,dname):
        return self.h5file.get(dname)
