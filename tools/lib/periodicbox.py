import h5
import flash
import flexi

probes = list()

def probe_flash(h5file):
    if 'sim info' in h5file.keys():
        names = h5file.get('sim info').dtype.names
        if 'flash version' in (x.lower() for x in names):
            return flash.File
    return None
probes.append(probe_flash)

# special handler for FLEXI which assumes special case
# of periodic box
def FlexiPeriodicBox(flexfp, meshfp=None, mode='r'):
    import os
    import hopr
    import pathlib as pl

    flexfp = pl.Path(flexfp)

    if meshfp is None:
        with h5.File(flexfp, mode='r') as fh:
            meshfp = flexfp.parent / pl.Path(fh.attrs['MeshFile'][0].decode('utf-8'))
    else:
        meshfp = pathlib.Path(meshfp)

    return flexi.File(flexfp, hopr.CartesianMeshFile(meshfp), mode)

def probe_flexi(h5file):
    if 'Program' in h5file.attrs.keys():
        pname = h5file.attrs['Program'][0].decode()
        if 'flexi' in pname.lower():
            return FlexiPeriodicBox
    return None
probes.append(probe_flexi)

def apply_probes(h5file):
    for probe in probes:
        result = probe(h5file)
        if result is not None:
            return result

def File(fpath, mode='r'):
    """Universal routine to probe the file type via introspection
    and calling the appropiate module on it."""

    with h5.File(fpath, mode='r') as fh:
        handler = apply_probes(fh) 

    if handler:
        return handler(fpath, mode=mode)
    
    raise NotImplementedError("No suitable handler for file '%s' found!" % fpath)
