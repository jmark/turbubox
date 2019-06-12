#!/usr/bin/env pyturbubox

import h5 as h5
import ulz as ulz
import numpy as np
import interpolate as itpl

class File(h5.File):
    def __init__(self, fpath, mode='r', **kwargs):
        super().__init__(fpath, mode)

        self.meta = self._get_metadata()

        for k in self.meta:
            setattr(self, k, self.meta[k])

    def _get_metadata(self):
        def _get_str(dname):
            dset = self.get(dname)
            return tuple(x.strip().decode() for x in dset[()])

        def _get_num(dname):
            dset = self.get(dname)
            return tuple(x for x in dset[()])

        keys = ()
        keys += _get_str('/meta/string/keys')
        keys += _get_str('/meta/integer/keys')
        keys += _get_str('/meta/float/keys')

        vals = ()
        vals += _get_str('/meta/string/values')
        vals += _get_num('/meta/integer/values')
        vals += _get_num('/meta/float/values')

        return dict(zip(keys,vals))

    @staticmethod
    def stitch(data):
        carpet = None

        for col in data:
            ribbon = None
            for row in col:
                if ribbon is None:
                    ribbon = row
                else:
                    ribbon = np.concatenate((ribbon,row),axis=0)
            if carpet is None:
                carpet = ribbon
            else:
                carpet = np.concatenate((carpet,ribbon),axis=1)
                
        return carpet

    @staticmethod
    def interpolate(patch,method=None,shape=None):

        if shape is None:
            shape = (patch.shape[0]*patch.shape[2],patch.shape[1]*patch.shape[3])

        image = np.zeros(shape)
        itpl.cells_to_image_titanic_patch_2d(patch,image,method=method)

        return image

    def get_var(self,varname,method=None,shape=None):
        if varname in ('dens','momx','momy','ener'):
            dpath = '/data/states/'+varname

            if method is None:
                return self.stitch(np.transpose(self.get(dpath)[()],(0,1,3,2)))

            return self.interpolate(self.get(dpath)[()],method,shape)

        if varname == 'blend':
            df = self.get('/data/hydro/blend')
            return tuple(self.stitch(np.transpose(df[:,:,i,:,:],(0,1,3,2))) for i in range(df.shape[2]))

        raise KeyError('Unknown varname: {0}'.format(varname))

if __name__ == '__main__':
    import sys
    fp = sys.argv[1]
    fh = File(fp)

    print(fh.program)
    print(fh.time)

    dens = fh.get_var('dens')
    blend = fh.get_var('blend')

    print(blend[1].shape)

    #import numpy as np
    #import matplotlib
    #matplotlib.rcParams.update({'font.size': 9})
    #import matplotlib.pyplot as plt

    ##fig = plt.figure(1,figsize=(geo[0]/dpi, geo[1]/dpi), dpi=dpi)
    #fig  = plt.figure()

    #plt.title('density: t = {:.2f}'.format(fh.time))
    #plt.imshow(
    #    data,
    #    #extent = extent,
    #    origin='lower',
    #    interpolation = None,
    #    cmap = plt.get_cmap('cubehelix'),
    #)
    #plt.colorbar()
    #plt.show()
