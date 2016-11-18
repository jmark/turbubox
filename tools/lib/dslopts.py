'''
Damn Small Commandline Arguments/Options Parse and Retrieval Agent

This is a nice KISS python module for parsing commandline arguments without
defining switches.

It builds a nice formmatted usage page automatically.

Example:

* conventional:
    import dslopts
    
    hdl = dslopts.Handler()
    hdl.arg(name='flshfile', desc='flash file path', type=Path, check=path_exists)
    args = hdl.parse() # returns a dictionary

    print(args['flshfile'])

* with context manager:
    import dslopts

    appendix = """  Methods are:

        0 -> without neighboring cells: n-th order interpolation (old version)
        1 -> with neighboring cells: (n+2)-th order interpolation (old version)   
        2 -> like '0' (new version)
        3 -> like '1' (new version)
        4 -> with neighboring cells which get averaged with boundary cells: n-th order interpolation"""

    with dslopts.Handler(scope=globals(),appendix=appendix) as hdl:
        hdl.arg(name='flshfile' ,desc='flash file path' ,type=Path  ,check=path_exists)
        hdl.arg(name='meshfile' ,desc='mesh file path'  ,type=Path  ,check=path_exists)
        hdl.arg(name='flexfile' ,desc='flexi file path' ,type=Path  ,check=path_exists)
        hdl.opt(name='method'   ,desc='method nr: 0-4'  ,type=int   ,check=meth_exists, default=3)

    print(meshfile)

    method == 2
'''

import sys

class Handler:
    def __init__(self, argsdict=None, scope=None, appendix=''):
        """
            argsdict    -> fill given dictionary with parsed parameters
            scope       -> install arguments as variables in given scope
            appendix    -> add appendix to the end of the usage message
        """

        self.ahandlers = list()
        self.ohandlers = list()
        self.arguments = {'_progname_': sys.argv[0]}
        self.helpkws   = 'help usage what how ?'.split()
        self.argsdict  = argsdict
        self.scope     = scope
        self.appendix  = appendix

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.parse()
        if self.argsdict:
            self.argsdict.update(self.arguments)
        if self.scope:
            self.install_into_scope(self.scope)

    def arg(self, name, desc='--', type=str, default=None, check=None):
        self.ahandlers.append({
            'name':         name,
            'desc':         desc,
            'type':         type,
            'check':        check,
            'default':      default,
        })

    def opt(self, name, desc='--', type=str, default=None, check=None):
        self.ohandlers.append({
            'name':         name,
            'desc':         desc,
            'type':         type,
            'check':        check,
            'default':      default,
        })

    def parse(self, ARGV=None):
        if not ARGV:
            ARGV = sys.argv[1:]

        ahdls = self.ahandlers # shortcut
        ohdls = self.ohandlers # shortcut

        # discriminate between arguments, options and ignored args
        argv = []
        optv = []
        ignv = []

        _tmp = argv
        for arg in ARGV:
            if arg == ':':
                _tmp = optv
                continue
            if arg == '::':
                _tmp = ignv
                continue
            _tmp.append(arg)

        # scan of help/usage arguments
        if any(x.lower() in self.helpkws for x in argv + optv):
            self.print_usage()
            sys.exit(1)

        # when arguments and associated handlers not match somthin' is fishy
        if len(argv) != len(ahdls):
            raise AssertionError(
                "Defined and given arguments list do not match up!\n\n"
                + "parameter list: '" + "' '".join(ARGV) + "'\n\n" + self.usage())

        # preset optional args
        for hdl in ohdls:
            self.arguments[hdl['name']] = hdl['default']

        def _parse(_argv, _hdls):
            for i, (arg, hdl) in enumerate(zip(_argv,_hdls),1):
                if arg is '-':
                    value = hdl['default']
                else:
                    # poor man's type checking
                    try:
                        value = hdl['type'](arg)
                    except ValueError:
                        raise ValueError(
                            "Parameter %d must be of type '%s'.\n\n" % (i, hdl['type'].__name__)
                            + "received:       '" + arg + "'\n"
                            + "all parameters: '" + "', '".join(ARGV) + "'\n\n" + self.usage())

                    # hook up user defined checking
                    if hdl['check']:
                        try:
                            value = hdl['check'](value)
                        except Exception:
                            raise AssertionError(
                                "Checking routine raised an error for parameter: %d\n\n" % (i)
                                + "received:       '" + arg + "'\n"
                                + "all parameters: '" + "', '".join(ARGV) + "'\n\n" + self.usage())

                self.arguments[hdl['name']] = value

        _parse(argv, ahdls)
        _parse(optv, ohdls)

        self.arguments['_ignored_'] = ignv
        return self.arguments

    def install_into_scope(self, scope):
        scope.update(self.arguments)

    def usage(self):
        ahdls = self.ahandlers # shortcut
        ohdls = self.ohandlers # shortcut
        hdls  = ahdls + ohdls

        hdName = 'name'
        hdType = 'type'
        hdDeft = 'default value'
        hdDesc = 'description'

        lenName = max([len(hdName)]+[len(x['name']) for x in hdls])
        lenType = max([len(hdType)]+[len(x['type'].__name__) for x in hdls])
        lenDeft = max([len(hdDeft)]+[len(str(x['default'])) for x in hdls])
        lenDesc = max([len(hdDesc)]+[len(x['desc']) for x in hdls])

        primer = "usage: %s [1] [2] ... : [n+1] [n+2] ... (optional args) :: ... (ignored args)\n\n" % self.arguments['_progname_']
        primer += "  * A hyphen '-' as argument activates default value.\n"
        primer += "  * Either '%s' triggers this help message." % "', '".join(self.helpkws)
        primer += " For more\n    information try: pydoc dslopts\n"
        primer += "\n"

        aheader = "   argn  | %-*s  | %-*s  | %-*s  | %-*s\n" % \
                    (lenName, hdName, lenType, hdType, lenDeft, hdDeft, lenDesc, hdDesc)
        stroke  = '  ' + '-' * (len(aheader)-2) + "\n"

        if ohdls:
            atable  = aheader + stroke 
            for i, hdl in enumerate(ahdls,1):
                atable += "    %3d  | %-*s  | %-*s  | %-*s  | %-*s\n" % (
                    i, lenName, hdl['name'], lenType, hdl['type'].__name__, lenDeft, str(hdl['default']), lenDesc, hdl['desc'])
        else:
            atable = "  No arguments defined.\n"

        if ohdls:
            otable = "\n   optn\n" + stroke
            for i, hdl in enumerate(ohdls,len(ahdls)+1):
                otable += "    %3d  | %-*s  | %-*s  | %-*s  | %-*s\n" % (
                    i, lenName, hdl['name'], lenType, hdl['type'].__name__, lenDeft, str(hdl['default']), lenDesc, hdl['desc'])
        else:
            otable = ""

        if self.appendix: otable += "\n"
        return primer + atable + otable + self.appendix

    def print_usage(self):
        print(self.usage(), file=sys.stderr)
