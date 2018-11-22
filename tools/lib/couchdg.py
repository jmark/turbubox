#!/usr/bin/env pyturbubox

import h5 as h5
import ulz as ulz

class File(h5.File):
    def __init__(self, fpath, mode='r', **kwargs):
        super().__init__(fpath, mode)

        self.meta = self._get_metadata()

        for k,v in self.meta.items():
            setattr(self, k, v)

        setattr(self, 'foo bar', 42)

    def _get_metadata(self):
        def _get_txt(dname):
            dset = self.get(dname)
            vals = tuple(x.strip().decode() for x in dset)
            keys = tuple(x.strip().decode() for x in dset.attrs.get('keys'))
            return tuple(filter(lambda kv: bool(kv[0]), zip(keys,vals)))

        def _get_num(dname):
            dset = self.get(dname)
            vals = dset
            keys = tuple(x.strip().decode() for x in dset.attrs.get('keys'))
            
            return tuple(filter(lambda kv: bool(kv[0]), zip(keys,vals)))

        return dict(_get_txt('meta_txt') + _get_num('meta_int') + _get_num('meta_flt'))

if __name__ == '__main__':
    import sys
    fp = sys.argv[1]
    fh = File(fp)

    print(fh.program)
    print(fh.time)

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
