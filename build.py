#!/usr/bin/env python


version_file = "QUICKFLASH_VERSION"
library_name = "QuickFlash"

library_file_name = "libquickflash"

machine_configfile = "machines.cfg"
settings_makefile = "settings.make"

lib_install_dir_default = "."
python_install_dir_default = "python"

make_cmd_default = "gmake"


# QuickFlash C++ files

quickflash_lib_dir = "src"
quickflash_swig_dir = "swig"


def read_version_file(filename) :

    ret_val = None

    ## Use the first non-empty, non-commented line

    file_ok = True

    f = None

    try :
        f = open(filename, "r")
    except :
        file_ok = False

    if (file_ok) :

        at_eof = False

        while ((not at_eof) and (ret_val is None)) :

            line = f.readline()

            if (line == "") :

                at_eof = True

            else :

                ## Remove both shell and C++-style single-line comments

                clean_line = line.split("#")[0].split("//")[0].strip()

                if (clean_line != "") :
                    ret_val = clean_line

    return ret_val

                
def try_cmd(command_str, error_message="Failed") :

    from sys import exit, stderr
    from os import system

    if (system(command_str) != 0) :
        stderr.write("\nError executing command [ " + command_str + " ] : "
                     + error_message + "\n\n")
        exit(-1)
    

def clean_up(make_cmd=make_cmd_default, use_python=False, python_only=False) :

    from sys import stderr
    from os import system

    if (not python_only) :
        cmd = "cd " + quickflash_lib_dir + " ; " + make_cmd + " clean"
        try_cmd(cmd, "Unable to clean library")

    if (use_python or python_only) :
        cmd = "cd " + quickflash_swig_dir + " ; " + make_cmd + " clean"
        try_cmd(cmd, "Unable to clean Python (SWIG) library")


def build_lib(make_cmd=make_cmd_default,
              lib_install_prefix=lib_install_dir_default,
              use_python=False, python_only=False,
              python_install_prefix=python_install_dir_default,
              static_only=False) :

    from os import makedirs
    from sys import stderr
    
    stderr.write("\nUsing make command " + make_cmd + "\n\n")
            
    install_path = lib_install_prefix

    if (install_path == "") :
        intall_path = "."

    include_path = install_path + "/include"
    lib_path = install_path + "/bin"

    python_path = python_install_prefix

    if (python_path == "") :
        python_path = "./python"

    python_path = python_path + "/quickflash"

    if (not python_only) :
        
        try :      
            makedirs(include_path)
            makedirs(lib_path)
        except OSError :
            cmd = "rm -f -r " + include_path + "/* " + lib_path + "/*"
            try_cmd(cmd, "Unable to remove include and lib contents")

        library_make_cmd = make_cmd

        if (static_only) :
            library_make_cmd += " " + library_file_name + ".a"

        cmd = "cd src ; " + library_make_cmd
        try_cmd(cmd, "Unable to build library")

        cmd = "cp src/*.hpp src/*.tcc " + include_path + "/."
        try_cmd(cmd, "Unable to copy library header files")

        cmd = "cp src/" + library_file_name + ".* " + lib_path + "/."

        try_cmd(cmd, "Unable to copy library object files")

    if (use_python or python_only) :

        try :
            makedirs(python_path)
        except OSError :
            cmd = "rm -f -r " + python_path + "/*"
            try_cmd(cmd, "Unable to remove python/QuickFlash contents")

        cmd = "cd swig ; " + make_cmd
        try_cmd(cmd, "Unable to build Python (SWIG) library")

        cmd = "cp -f swig/quickflash.py swig/_quickflash.so " \
              + "swig/python/* " + python_path + "/."
        try_cmd(cmd, "Unable to copy Python (SWIG) library")


def main() :
    
    from sys import argv, exit, stdout, stderr
    from optparse import OptionParser
    from config import read_config, MachineSetup, make_key, so_ext_key, \
         lib_install_key, python_install_key, read_makefile

    ## Get the version string

    version_string = read_version_file(version_file)

    if (version_string is None) :
        version_string = ""

    ## Set up option parser

    usage = "Usage: %prog [ options ]"
    version = "%prog for " + library_name + " " + version_string

    argparser = OptionParser(usage, version=version)

    machine_config_str = str("machine-config")
    argparser.add_option("--" + machine_config_str, type="string",
                         dest="machine_config",
                         help=("Machine configuration file "
                               + "(default: " + machine_configfile + ")"))
    argparser.set_defaults(machine_config=machine_configfile)

    machine_context_str = str("machine")
    argparser.add_option("--" + machine_context_str, type="string",
                         dest="machine_context",
                         help=("Generates new settings.make file for given "
                               + "machine context "
                               + "(default: use existing settings.make)"))
    argparser.set_defaults(machine_context_str=None)

    machine_help_str = str("machine-help")
    argparser.add_option("--" + machine_help_str, action="store_true",
                         dest="machine_help",
                         help="List available machine descriptions")
    argparser.set_defaults(machine_help=False)

    screen_only_str = str("screen-only")
    argparser.add_option("--" + screen_only_str, action="store_true",
                         dest="screen_only",
                         help="Send settings.make contents to screen only")
    argparser.set_defaults(screen_only=False)

    enable_python_str = str("enable-python")
    argparser.add_option("--" + enable_python_str, action="store_true",
                         dest="enable_python",
                         help="Build Python bindings")
    argparser.set_defaults(enable_python=False)

    python_only_str = str("python-only")
    argparser.add_option("--" + python_only_str, action="store_true",
                         dest="python_only",
                         help=("Build or clean only Python bindings "
                               + "(assumes C++ library already built)"))
    argparser.set_defaults(python_only=False)

    clean_str = str("clean")
    argparser.add_option("--" + clean_str, action="store_true",
                         dest="clean_build",
                         help="Clean up build directores")
    argparser.set_defaults(clean_build=False)

    clean_python_str = str("clean-python")
    argparser.add_option("--" + clean_python_str, action="store_true",
                         dest="clean_python",
                         help="Shortcut for --python-only --clean")
    argparser.set_defaults(clean_python=False)

    clean_first_str = str("clean-first")
    argparser.add_option("--" + clean_first_str, action="store_true",
                         dest="clean_first",
                         help="Run clean before other operations")
    argparser.set_defaults(clean_python=False)

    makefile_only_str = str("makefile-only")
    argparser.add_option("--" + makefile_only_str, action="store_true",
                         dest="makefile_only",
                         help="Generate settings.make only")
    argparser.set_defaults(makefile_only=False)

    make_cmd_str = str("make-cmd")
    argparser.add_option("--" + make_cmd_str, type="string",
                         dest="make_cmd",
                         help=("Specify make command "
                               + "(overrides machine settings)"))
    argparser.set_defaults(make_cmd=None)

    static_only_str = str("static-only")
    argparser.add_option("--" + static_only_str, action="store_true",
                         dest="static_only",
                         help=("Do not generate shared object version of "
                               + library_file_name
                               + " (does not affect Python bindings)"))
    argparser.set_defaults(static_only=False)


    ## Set up the control variables

    config_filename = None
    machine_context = None

    make_cmd = None

    static_only = None

    build_library = True  ## Default behavior
    clean_library = False

    no_clean_build = False  ## If set, do not allow builds or cleans

    machine_help = False

    show_makefile = False
    write_makefile = False

    use_python = False
    python_only = False


    ## Process the command line input
    
    opts, args = argparser.parse_args()

    config_filename = opts.machine_config
    machine_context = opts.machine_context

    make_cmd = opts.make_cmd

    static_only = opts.static_only

    if (opts.machine_help) :
        machine_help = True
        build_library = False
        no_clean_build = True
    
    if (opts.screen_only) :
        show_makefile = True
        build_library = False
        no_clean_build = True

    if (opts.makefile_only) :
        makefile_only = True
        build_library = False
        no_clean_build = True

    if (opts.enable_python) :
        use_python = True

    if (opts.python_only) :
        use_python = True
        python_only = True

    if (opts.clean_build) :

        if (no_clean_build) :
            stderr.write("\nError: Options incompatible with clean operation"
                         + "\n\n")
            exit(-1)

        build_library = False
        clean_library = True
        
    if (opts.clean_python) :

        if (no_clean_build) :
            stderr.write("\nError: Options incompatible with clean operation"
                         + "\n\n")
            exit(-1)

        python_only = True
        build_library = False
        clean_library = True

    if (opts.clean_first) :

        if (no_clean_build) :
            stderr.write("\nError: Options incompatible with clean operation"
                         + "\n\n")
            exit(-1)
        
        clean_library = True


    ## Exit if unrecognized arguments exist

    if (len(args) != 0) :
        
        err = "Unrecognized arguments " + str(args)

        stderr.write("\nError: " + err + "\n\n")
        exit(-1)


    # Read the context list

    mach_setup = MachineSetup()

    lib_path = str(lib_install_dir_default)
    python_path = str(python_install_dir_default)

    if (machine_help) :
        
        context_list = read_config(config_filename)

        context_names = context_list.keys()
        context_names.sort()

        stderr.write("\nMachine context names\n\n")

        for name in context_names :
            if (name != "") :

                stderr.write("  " + name + "\n")
                
                help_str = context_list[name].get_help_string()

                if (help_str is not None) :
                    stderr.write("      " + help_str + "\n")

                stderr.write("\n")

        stderr.write("\n")
        
    elif (machine_context is not None) :

        # Read the setup

        context_list = read_config(config_filename)
    
        if (context_list.has_key(machine_context)) :

            mach_setup.reset(context_list[machine_context])
        
        else :
        
            err = "Machine context [ " + machine_context \
                  + " ] not present in file [ " + config_filename + " ]"
        
            stderr.write("\nError: " + err + "\n\n")
            exit(-1)

        if (make_cmd is None) :
            
            make_cmd = str(make_cmd_default)

            if (mach_setup.has_key(make_key)) :
                make_cmd = str(mach_setup[make_key])

        if (mach_setup.has_key(lib_install_key)) :
            lib_path = mach_setup[lib_install_key]

        if (mach_setup.has_key(python_install_key)) :
            python_path = mach_setup[python_install_key]

        # Build settings.make
   
        if (show_makefile) :

            # Send the makefile contents to stdout

            stdout.write("\n")
        
            mach_setup.write_settings_makefile(context_name=machine_context)
        
        else :

            ## Write the makefile
        
            mach_setup.write_settings_makefile(settings_makefile,
                                               context_name=machine_context)

    elif (show_makefile) :

        ## Read the settings makefile

        stderr.write("\nWarning: No machine context specified -- using "
                     + settings_makefile + "\n\n")

        settings_config = read_makefile(settings_makefile)

        mach_setup.reset(settings_config)

        mach_setup.write_settings_makefile()


    elif (make_cmd is None) :
        
        ## See if the settings makefile has a value for make command

        settings_config = read_makefile(settings_makefile)

        if (settings_config.has_key(make_key)) :

            make_cmd_words = settings_config[make_key]

            num_words = len(make_cmd_words)

            if (num_words < 1) :
                stderr.write("\nError: Empty make command string")
                exit(-1)

            make_cmd = str(make_cmd_words[0])

            if (num_words > 1) :
                for index in range(1, num_words) :
                    make_cmd += " " + str(make_cmd_words[index])

        else :
            
            make_cmd = str(make_cmd_default)


    ## Clean and/or build the library

    if (not no_clean_build) :
            
        if (clean_library) :
            clean_up(make_cmd, use_python, python_only)
            
        if (build_library) :
            build_lib(make_cmd, lib_path, use_python, python_only, python_path,
                      static_only)


if (__name__ == "__main__") :
    main()
