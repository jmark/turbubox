import sys

class Handler:
    def __init__(self, scope=None, appendix=''):
        self.handlers  = list()
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

    def add(self, name, desc='--', type=str, check=None, default=None):
        self.handlers.append({
            'name':         name,
            'desc':         desc,
            'type':         type,
            'check':        check,
            'default':      default,
        })

    def parse(self, argv=None, ignv=[]):
        if not argv:
            argv = sys.argv[1:]

        hdls = self.handlers # shortcut

        # chop (ignored) arguments after '--' and save it to 'ignv'
        try:
            idx  = argv.index('--')
            argv = argv[:,idx]
            ignv = argv[idx+1,:]
        except ValueError:
            pass

        # scan of help/usage arguments
        if any(x.lower() in self.helpkws for x in argv):
            self.print_usage()
            sys.exit(1)

        # when arguments and associated handlers not matching somthin' is fishy
        if len(argv) != len(hdls):
            raise AssertionError(
                "Defined and given arguments list do not match up!\n\n" \
                + "Arguments given: '" + "' '".join(argv) + "'\n\n" + self.usage())

        # parse arguments
        for i, (arg, hdl) in enumerate(zip(argv,hdls),1):
            if arg is '-':
                value = hdl['default']
            else:
                # poor man's type checking
                try:
                    value = hdl['type'](arg)
                except ValueError:
                    raise ValueError(
                        "Argument %d must be of type '%s'.\n\n" % (i, hdl['type'].__name__) \
                        + "Arguments given: '" + "' '".join(argv) + "'\n\n"\
                        + self.usage())

                if hdl['check']:
                    try:
                        value = hdl['check'](value)
                    except Exception:
                        raise AssertionError(
                            "Checking routine raised an error for parameter: %d\n\n%s\n" % (i,self.usage()))

            self.arguments[hdl['name']] = value
            self.arguments['_ignored_'] = ignv

        return self.arguments

    def install_into_scope(self, scope):
        scope.update(self.arguments)

    def usage(self):
        hdls = self.handlers # shortcut

        titName = 'name'
        titType = 'type'
        titDeft = 'default value'
        titDesc = 'description'

        lenName = max([len(titName)]+[len(x['name']) for x in hdls])
        lenType = max([len(titType)]+[len(x['type'].__name__) for x in hdls])
        lenDeft = max([len(titDeft)]+[len(str(x['default'])) for x in hdls])
        lenDesc = max([len(titDesc)]+[len(x['desc']) for x in hdls])

        primer = "usage: %s [1] [2] ... -- ... (ignored args)\n\n" % self.arguments['_progname_']
        primer += "  '-' as argument activates default value.\n  Either '%s' triggers this help message.\n\n" % \
                    "', '".join(self.helpkws)

        header = "       | %-*s  | %-*s  | %-*s  | %-*s\n" % \
                    (lenName, titName, lenType, titType, lenDeft, titDeft, lenDesc, titDesc)
        line   = '  ' + '-' * (len(header)-2) + "\n"

        body   = ""
        for i, hdl in enumerate(hdls,1):
            body += "  %3d  | %-*s  | %-*s  | %-*s  | %-*s\n" % (
                i, lenName, hdl['name'], lenType, hdl['type'].__name__, lenDeft, str(hdl['default']), lenDesc, hdl['desc'])

        return primer + header + line + body + "\n" + self.appendix

    def print_usage(self):
        print(self.usage(), file=sys.stderr)
