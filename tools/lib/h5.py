import h5py
import sys

class H5File(h5py.File):
    def __init__(self,fpath,mode='r'):
        try:
            super().__init__(fpath,mode)
        except Exception as e:
            print("%s" % e, file=sys.stderr)
            sys.exit(1)
    
    def get(self,dname):
        db = super().get(dname)
        if db:
            return db
        else:
            raise KeyError(dname)

File = H5File

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
