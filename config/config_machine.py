# Configuration key names

make_key = "MAKE"

use_package_key = "USE_PACKAGE"

lib_install_key = "LIB_INSTALL_DIR"
python_install_key = "PYTHON_INSTALL_DIR"

cxx_key = "CXX"
ld_key = "LD"

ccflags_key = "CCFLAGS"
ldflags_key = "LDFLAGS"

ccflags_shared_key = "CCFLAGS_SHARED"
ldflags_shared_key = "LDFLAGS_SHARED"

ccflags_noshared_key = "CCFLAGS_NOSHARED"
ldflags_noshared_key = "LDFLAGS_NOSHARED"

ccflags_hdf5_key = "CCFLAGS_HDF5"
ldflags_hdf5_key = "LDFLAGS_HDF5"

so_ext_key = "SO_EXT"

ranlib_key = "RANLIB"
ar_key = "AR"
arflags_key = "ARFLAGS"

swig_key = "SWIG"

cxx_python_key = "CXX_PYTHON"
ld_python_key = "LD_PYTHON"

ccflags_python_key = "CCFLAGS_PYTHON"
ldflags_python_key = "LDFLAGS_PYTHON"

ccflags_shared_python_key = "CCFLAGS_SHARED_PYTHON"
ldflags_shared_python_key = "LDFLAGS_SHARED_PYTHON"

ccflags_gsl_key = "CCFLAGS_GSL"
ldflags_gsl_key = "LDFLAGS_GSL"

ccflags_magick_key = "CCFLAGS_MAGICK"
ldflags_magick_key = "LDFLAGS_MAGICK"


# Default configuration

make_default = "gmake"

use_package_default = ""

lib_install_default = "."
python_install_default = "./python"

cxx_default = "CC"
ld_default = None  # If not specified, use CXX

ccflags_default = "-g -O2"
ldflags_default = "-g"

ccflags_shared_default = ""
ldflags_shared_default = "-shared"

ccflags_noshared_default = ""
ldflags_noshared_default = ""

ccflags_hdf5_default = "-I /usr/include"
ldflags_hdf5_default = "-L/usr/lib -lhdf5"

so_ext_default = "so"

ranlib_default = "ranlib"
ar_default = "ar"
arflags_default = "rsc"

swig_default = "swig"

cxx_python_default = None  # If not specified, use CXX
ld_python_default = None  # If not specified, use CXX_PYTHON

ccflags_shared_python_default = None  # If not specified, use CCFLAGS_SHARED
ldflags_shared_python_default = None  # If not specified, use LDFLAGS_SHARED

ccflags_python_default = "-I /usr/include/python"
ldflags_python_default = ""

ccflags_gsl_default = "-I /usr/include"
ldflags_gsl_default = "-L/usr/lib -lgsl -lgslcblas"

ccflags_magick_default = ""
ldflags_magick_default = ""


class MachineSetup :

    def __init__(self, context_name=None, help_string=None) :

        self.__params = None

        self.__context_name = context_name
        self.__help_str = help_string

        self.__reset_default()


    def __reset_default(self) :

        self.__params = dict()

        self.__params[make_key] = make_default

        self.__params[use_package_key] = use_package_default

        self.__params[lib_install_key] = lib_install_default
        self.__params[python_install_key] = python_install_default

        self.__params[cxx_key] = cxx_default
        self.__params[ld_key] = ld_default

        self.__params[ccflags_key] = ccflags_default
        self.__params[ldflags_key] = ldflags_default

        self.__params[ccflags_shared_key] = ccflags_shared_default
        self.__params[ldflags_shared_key] = ldflags_shared_default

        self.__params[ccflags_noshared_key] = ccflags_noshared_default
        self.__params[ldflags_noshared_key] = ldflags_noshared_default

        self.__params[ccflags_hdf5_key] = ccflags_hdf5_default
        self.__params[ldflags_hdf5_key] = ldflags_hdf5_default

        self.__params[so_ext_key] = so_ext_default

        self.__params[ranlib_key] = ranlib_default
        self.__params[ar_key] = ar_default
        self.__params[arflags_key] = arflags_default

        self.__params[swig_key] = swig_default

        self.__params[cxx_python_key] = cxx_python_default
        self.__params[ld_python_key] = ld_python_default

        self.__params[ccflags_shared_python_key] \
                                  = ccflags_shared_python_default
        self.__params[ldflags_shared_python_key] \
                                  = ldflags_shared_python_default

        self.__params[ccflags_python_key] = ccflags_python_default
        self.__params[ldflags_python_key] = ldflags_python_default

        self.__params[ccflags_gsl_key] = ccflags_gsl_default
        self.__params[ldflags_gsl_key] = ldflags_gsl_default

        self.__params[ccflags_magick_key] = ccflags_magick_default
        self.__params[ldflags_magick_key] = ldflags_magick_default


    def reset(self, config_context) :

        self.__reset_default()

        self.__context_name = config_context.get_context_name()
        self.__help_str = config_context.get_help_string()

        for key in config_context.get_keys() :
            self.__params[key] = config_context.get_value_str(key)


    def get_context_name() :
        return self.__context_name
    

    def get_help_string() :
        return self.__help_str
    

    def has_key(self, key) :
        return self.__params.has_key(key)
    

    def get_param(self, key) :

        ret_val = str("")

        if (self.__params.has_key(key)) :
            ret_val = str(self.__params[key])

        return ret_val


    def set_param(self, key, value) :
        self.__params[key] = str(value)


    def __getitem__(self, key) :
        return self.get_param(key)


    def __setitem__(self, key, value) :
        self.set_param(key, value)


    def __check_compilers(self) :

        # Check C++ compiler

        if (self.__params[cxx_key] is None) :
            raise ValueError("No C++ compiler specified")

        if (self.__params[cxx_key] == "") :
            raise ValueError("Empty C++ compiler string")

        # Check the linker

        if (self.__params[ld_key] is None) :
            self.__params[ld_key] = str(self.__params[cxx_key])
        elif (self.__params[ld_key] == "") :
            raise ValueError("Empty linker command string")

        # Check the Python C++ compiler

        if (self.__params[cxx_python_key] is None) :
            self.__params[cxx_python_key] = str(self.__params[cxx_key])
        elif (self.__params[cxx_python_key] == "") :
            raise ValueError("Empty C++ compiler string for Python")

        # Check the Python linker

        if (self.__params[ld_python_key] is None) :
            self.__params[ld_python_key] = str(self.__params[cxx_python_key])
        elif (self.__params[ld_python_key] == "") :
            raise ValueError("Empty linker command string for Python")

        # Check the Python shared flags

        if (self.__params[ccflags_shared_python_key] is None) :
            self.__params[ccflags_shared_python_key] \
                                      = str(self.__params[ccflags_shared_key])

        if (self.__params[ldflags_shared_python_key] is None) :
            self.__params[ldflags_shared_python_key] \
                                      = str(self.__params[ldflags_shared_key])


    def write_config_context(self, filename=None, context_name=None,
                             help_string=None) :

        from sys import stdout
        from config_file import help_string_marker, left_bracket_marker, \
             right_bracket_marker

        self.__check_compilers()

        # SAVE PREVIOUS VERSION????

        f = stdout

        if (filename is not None) :
            f = file(filename, "a")

        f.write("\n")

        context_title = context_name

        if (context_title is None) :
            context_title = self.__context_name

        help_str = help_string

        if (help_str is None) :
            help_str = self.__help_str
        
        if (context_title is not None) :

            title_string = context_title

            if (help_str is not None) :
                title_string += " " + help_string_marker + " " + help_str
                
            f.write(left_bracket_marker + " " + str(title_string) + " "
                    + right_bracket_marker + "\n")

        for key in self.__params :

            param_str = self.__params[key]

            if (param_str is None) :
                param_str = ""
 
            f.write(str(key) + " " + str(param_str) + "\n")


    def write_settings_makefile(self, filename=None, context_name=None,
                                help_string=None) :

        from sys import stdout

        self.__check_compilers()

        # SAVE PREVIOUS VERSION???

        f = stdout

        if (filename is not None) :
            f = file(filename, "w")

        context_title = context_name

        if (context_title is None) :
            context_title = self.__context_name

        help_str = help_string

        if (help_str is None) :
            help_str = self.__help_str

        header_present = False
        
        if (context_title is not None) :
            f.write("# Makefile for setup [ " + str(context_title) + " ]\n")
            header_present = True

        if (help_str is not None) :
            f.write("# " + str(help_str) + "\n")
            header_present = True

        if (header_present) :
            f.write("\n")

        for key in self.__params :
            
            param_str = self.__params[key]

            if (param_str is None) :
                param_str = ""
            
            ## The Python on Unclassified Purple somehow redefines the empty
            ## string("") as a PyCObject.  Let's use a non-empty string for
            ## now.
            ##
            ## NOTE: The definition change occurs during the
            ## f = file(filename, "w") statement above

            if (param_str == "") :
                param_str = " "

            ## Add the external package usage flags

            if (key == ccflags_key) :

                use_param = self.__params[use_package_key]

                ## Keep IBM Python from complaining ...

                if (use_param != "") :

                    use_words = use_param.split()
                    
                    for word in use_words :

                        param_str += " -DUSE_" + word
                
            f.write(str(key) + " = " + str(param_str) + "\n")
