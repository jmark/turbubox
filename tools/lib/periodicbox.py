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
def FlexiPeriodicBox(srcfp, meshfp=None, mode='r'):
    import os
    import hopr
    if meshfp is None:
        h5file = h5.File(srcfp, 'r')
        meshfp = h5file.attrs['MeshFile'][0].decode('utf-8')
        h5file.close()
        oldcwd = os.path.realpath(os.curdir)
        #print(os.path.dirname(os.path.realpath(srcfp)))
        os.chdir(os.path.dirname(os.path.realpath(srcfp)))
        meshfp = os.path.realpath(meshfp)
        os.chdir(oldcwd)

    return flexi.File(srcfp, hopr.CartesianMeshFile(meshfp), mode)

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

    h5file = h5.File(fpath, mode='r')
    handler = apply_probes(h5file) 
    h5file.close()

    if handler:
        return handler(fpath, mode=mode)
    
    raise NotImplementedError("No suitable handler for file '%s' found!" % fpath)
