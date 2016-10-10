import h5py as h5
import sys

class H5File(h5.File):
    def __init__(self,fpath,mode):
        try:
            super().__init__(fpath,mode)
        except Exception as e:
            print("%s" % e, file=sys.stderr)
            sys.exit(1)
    
    def get(self,dname):
        db = super(H5File,self).get(dname)
        if db:
            return db
        else:
            raise KeyError(dname)

    @staticmethod
    def dataset_to_dict(dataset):
        """
        Transform a key/value dataset to a pythonic dictionary

        Parameters
        ----------
        dataset : str
            dataset name

        Returns
        -------
        dictionary
        """
        return dict([(k.strip().decode(),v) for (k,v) in dataset])
