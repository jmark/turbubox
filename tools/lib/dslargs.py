import sys

class Handler:
    def __init__(self, scope=None, appendix=''):
        self.ahandlers  = list()
        self.ohandlers = list()
        self.arguments = {'_progname_': sys.argv[0]}
        self.helpkws   = 'help usage what how ?'.split()
        self.scope     = scope
        self.appendix  = appendix

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.parse()
        if self.scope:
            self.install_into_scope(self.scope)

    def arg(self, name, desc='--', type=str, check=None, default=None):
        self.ahandlers.append({
            'name':         name,
            'desc':         desc,
            'type':         type,
            'check':        check,
            'default':      default,
        })

    def opt(self, name, desc='--', type=str, check=None, default=None):
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

        argv = []
        optv = []
        ignv = []

        ll = argv
        for arg in ARGV:
            if arg == ':':
                ll = optv
                continue
            if arg == '::':
                ll = ignv
                continue
            ll.append(arg)

        # scan of help/usage arguments
        if any(x.lower() in self.helpkws for x in argv + optv):
            self.print_usage()
            sys.exit(1)

        # when arguments and associated handlers not matching somthin' is fishy
        if len(argv) != len(ahdls):
            raise AssertionError(
                "Defined and given arguments list do not match up!\n\n"
                + "parameter list: '" + "' '".join(ARGV) + "'\n\n" + self.usage())

        # parse obligatory arguments
        for i, (arg, hdl) in enumerate(zip(argv,ahdls),1):
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

                if hdl['check']:
                    try:
                        value = hdl['check'](value)
                    except Exception:
                        raise AssertionError(
                            "Checking routine raised an error for parameter: %d\n\n" % (i)
                            + "received:       '" + arg + "'\n"
                            + "all parameters: '" + "', '".join(ARGV) + "'\n\n" + self.usage())


            self.arguments[hdl['name']] = value

        # preset optional args
        for hdl in ohdls:
            self.arguments[hdl['name']] = hdl['default']

        # parse optional arguments
        for i, (arg, hdl) in enumerate(zip(optv,ohdls),len(argv)+1):
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

                if hdl['check']:
                    try:
                        value = hdl['check'](value)
                    except Exception:
                        raise AssertionError(
                            "Checking routine raised an error for parameter: %d\n\n" % (i)
                            + "received:       '" + arg + "'\n"
                            + "all parameters: '" + "', '".join(ARGV) + "'\n\n" + self.usage())

            self.arguments[hdl['name']] = value

        self.arguments['_ignored_'] = ignv
        return self.arguments

    def install_into_scope(self, scope):
        scope.update(self.arguments)

    def usage(self):
        ahdls = self.ahandlers # shortcut
        ohdls = self.ohandlers # shortcut
        hdls  = ahdls + ohdls

        titName = 'name'
        titType = 'type'
        titDeft = 'default value'
        titDesc = 'description'

        lenName = max([len(titName)]+[len(x['name']) for x in hdls])
        lenType = max([len(titType)]+[len(x['type'].__name__) for x in hdls])
        lenDeft = max([len(titDeft)]+[len(str(x['default'])) for x in hdls])
        lenDesc = max([len(titDesc)]+[len(x['desc']) for x in hdls])

        primer = "usage: %s [1] [2] ... : [n+1] [n+2] ... (optional args) :: ... (ignored args)\n\n" % self.arguments['_progname_']
        primer += "  * A hyphen '-' as argument activates default value.\n"
        primer += "  * Either '%s' triggers this help message.\n\n" % "', '".join(self.helpkws)

        aheader = "   argn  | %-*s  | %-*s  | %-*s  | %-*s\n" % \
                    (lenName, titName, lenType, titType, lenDeft, titDeft, lenDesc, titDesc)
        line    = '  ' + '-' * (len(aheader)-2) + "\n"

        atable  = aheader + line
        for i, hdl in enumerate(ahdls,1):
            atable += "    %3d  | %-*s  | %-*s  | %-*s  | %-*s\n" % (
                i, lenName, hdl['name'], lenType, hdl['type'].__name__, lenDeft, str(hdl['default']), lenDesc, hdl['desc'])

        if ohdls:
            otable = "\n   optn\n" + line
            for i, hdl in enumerate(ohdls,len(ahdls)+1):
                otable += "    %3d  | %-*s  | %-*s  | %-*s  | %-*s\n" % (
                    i, lenName, hdl['name'], lenType, hdl['type'].__name__, lenDeft, str(hdl['default']), lenDesc, hdl['desc'])
        else:
            otable = ""

        if self.appendix: otable += "\n"
        return primer + atable + otable + self.appendix

    def print_usage(self):
        print(self.usage(), file=sys.stderr)
