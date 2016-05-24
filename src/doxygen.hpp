// C++ header/Doxygen format file

/*
  By Nathan C. Hearn
     January 2, 2009

  Doxygen mainpage text for the QuickFlash library.
*/

/*!

\mainpage QuickFlash Documentation

Visit the <a href="http://quickflash.sourceforge.net">QuickFlash Homepage</a>


\author Nathan C. Hearn<br>Contributors: George C. Jordan, Dean Townsley, John Norris
\version 1.0.0


\section sec_intro Introduction

%QuickFlash is an open-source C++ analysis library for 
Paramesh-structured HDF5 data files produced by
<a href="http://flash.uchicago.edu">Flash</a>, an Eulerian
hydrodynamics code with adaptive mesh refinement developed by the
<a href="http://flash.uchicago.edu">ASC/Alliances Center for 
Astrophysical Thermonuclear Flashes</a> (the "Flash Center") at the 
<a href="http://www.uchicago.edu">University of Chicago</a>.  
%QuickFlash provides high-performance data access and processing 
routines for large, multi-dimensional datasets in environments where 
memory and access to other software may be limited, such as on 
desktop computers or parallel supercomputer nodes.

The core of the %QuickFlash library focuses on direct access to data
stored in AMR blocks, as well as efficient means of locating the 
desired data by spatial coordinates.  %QuickFlash uses caches and read-ahead
buffers to provide fast, on-demand loading of block data while limiting
the demands on memory and disk access.  The API includes simple access
functions to AMR tree information for the purpose of quickly locating
block indexes, as well as for traversing the tree by branch and node.

%QuickFlash includes higher level functionality, including slicing, 
binning, and isosurface construction routines.  When compiled with the optional
<a href="http://www.imagemagick.org">ImageMagick</a> C++ bindings (Magick++), 
%QuickFlash can be used to generate
images from data with continous, user-defined colormaps.  Example 
applications provided with the library include slicing, raytracing,
fractal boxcounting, and isosurface utilities.

The %QuickFlash library is released under 
<a href="http://www.gnu.org/licenses/gpl-2.0.html">Version 2</a> of the
<a href="http://www.gnu.org/licenses">GNU General Public License</a>.  
In short, this license gives you freedom to use, modify, and 
distribute the software as you wish, as long as any 
distributed versions of the code are offered under the same terms
with prominant notices of any modifications to the code.

If you find %QuickFlash useful, please register your name and email 
address with the maintainer according to the README file included 
with the library.  Your information will not be shared with third 
parties, but knowing the number of active users of the code will help 
in supporting continued development of the code.


\section sec_comp Compiling QuickFlash

The build system is a Python script (<tt>build.py</tt>) that performs the following operations:

<ol>
  <li>Generates a Makefile header (<tt>settings.make</tt>) that includes definitions for compiling the %QuickFlash library and related programs</li>

  <li>Builds the library in the <tt>src_lib</tt> directory</li>

  <li>(<em>Optional</em>) Builds the %QuickFlash Python bindings in the <tt>swig</tt> directory</li>

  <li>Copies the library header files and object files to the installation directories</li>
</ol>

Brief documentation about the various options for <tt>build.py</tt> can be found by issuing the command:

\verbatim
./build.py --help
\endverbatim

More details can be found in section \ref subsec_comp_build_py.


\subsection subsec_comp_config_build Configuring the Build System

By default, the <tt>build.py</tt> script uses information in the <tt>machines.cfg</tt> file to construct <tt>settings.make</tt>.  A setup for a given build includes specifications for the compiler and external libraries (such as HDF5).  The <tt>machines.cfg</tt> file can hold multiple setup specifications (called <em>machine contexts</em>), each denoted by a label.

A machine context is specified by starting with the label encased in brackets, as in

\verbatim
[ my_context ]
\endverbatim

Spaces between the first bracket and the start of label and those between the end of the label and the last bracket are ignored.  An optional description can be added by placing a colon followed by the text after the context name, such as

\verbatim
[ my_context : my description ]
\endverbatim

In general, spaces are not recommended for use within the context name string &mdash; underscores or dashes work well as substitues &mdash; but they are welcome in the description text.

Definitions follow the label, where each line begins with a keyword followed by the definition associated with that keyword. Mulitiple-line entries can be constructed by placing a "backslash" character (<tt>\\</tt>) at the end of each line that is to be extended.

The line

\verbatim
CXX g++
\endverbatim

assigns the value <tt>g++</tt> to the symbol <tt>CXX</tt>.  The keyword and value will be placed in the <tt>settings.make</tt> file using the same case.  The definition above will yield the line

\verbatim
CXX = g++
\endverbatim

in the <tt>settings.make</tt> file.  (NOTE: At this time, there is no guarantee that entries in a <tt>machines.cfg</tt> context will have their order preserved when <tt>settings.make</tt> is built.)  Comments can be added to <tt>machines.cfg</tt> by using the comment character <tt>#</tt>.  Any text between a <tt>#</tt> and the end of the line is ignored.  Blank lines are also ignored in <tt>machines.cfg</tt>.

Any desired definitions can be placed within the context, all of which will be transferred to <tt>settings.make</tt>.  (However, keywords may not contain spaces or control characters.)  In order to construct the %QuickFlash library, the following keywords must be defined:

<ul>
  <li>
    <tt>MAKE</tt> &mdash; The name of the executable used for running 
    makefiles (usually <em>gmake</em> or <em>make</em>)
  </li>
  <li><tt>CXX</tt> &mdash; The C++ compiler</li>
  <li><tt>LD</tt> &mdash; The linker (uses <tt>CXX</tt> setting by default)</li>
  <li><tt>CCFLAGS</tt> &mdash; Essential compilation flags</li>
  <li><tt>LDFLAGS</tt> &mdash; Essential linker flags</li>
  <li><tt>CCFLAGS_SHARED</tt> &mdash; Additional compilation flags when building shared libraries</li>
  <li><tt>LDFLAGS_SHARED</tt> &mdash; Additional linker flags when building shared libraries</li>
  <li><tt>CCFLAGS_NOSHARED</tt> &mdash; Additional compilation flags when not building shared libraries</li>
  <li><tt>LDFLAGS_NOSHARED</tt> &mdash; Additional linker flags when not building shared libraries</li>
  <li><tt>SO_EXT</tt> &mdash; Default system-specific extension used for shared libraries, without the leading period (usually <em>so</em> or <em>dylib</em>)</li>
  <li><tt>AR</tt> &mdash; Name of the executable used for building static archives from object files (usually <em>ar</em>)</li>
  <li><tt>ARFLAGS</tt> &mdash; Flags used when building static archives</li>
  <li><tt>RANLIB</tt> &mdash; Name of the executable used for adding indexes to static archives (usually <em>ranlib</em>)</li>
  <li><tt>CCFLAGS_HDF5</tt> &mdash; Additional compilation flags for the HDF5 library</li>
  <li><tt>LDFLAGS_HDF5</tt> &mdash; Additional linker flags for the HDF5 library</li>
  <li><tt>LIB_INSTALL_DIR</tt> &mdash; Path to directory where the library <tt>include</tt> and <tt>lib</tt> directories will be installed (NOTE: Custom installation directories not supported at this time).</li>
</ul>

Many of these keywords have default values; for any required keywords not found in the <tt>machines.cfg</tt> file, the default will be used if it exists.


\subsubsection subsubsec_build_python_settings Build Settings for Python Bindings

To build the Python bindings, the following additional keywords are needed:

<ul>
  <li><tt>SWIG</tt> &mdash; The name of the SWIG interface builder executable (usually <em>swig</em>); should be version 1.3.29 or later</li>
  <li><tt>CCFLAGS_PYTHON</tt> &mdash; Additional compiler flags when building Python packages (usually an include search path flag indicating the location of <tt>Python.h</tt>)</li>
  <li><tt>LDFLAGS_PYTHON</tt> &mdash; Additional compiler flags when building Python packages</li>
  <li><tt>PYTHON_INSTALL_DIR</tt> &mdash; Path to directory where the QuickFlash Python pakcage will be installed (NOTE: Not supported at this time)</li>
</ul>

Some of these keywords also have defaults.  Furthermore, the swig directory contains pre-built versions of the output from SWIG (<tt>%QuickFlash.py</tt> and <tt>python_quickflash_wrap.cxx</tt>), and make will use them if they are recent enough; thus, having the SWIG executable present may not be necessary to build the Python bindings for %QuickFlash.


\subsubsection subsubsec_build_opt_comps Optional Components

The %QuickFlash build system has support for optional code components
for both the library and end user programs.  The <tt>USE_PACKAGE</tt>
keyword in <tt>machines.cfg</tt> setups allows the user to specify
codes for one or more optional components.  Placing the code 
<tt>EXAMPLE_PKG</tt> in a setup via the line

\verbatim
USE_PACKAGE EXAMPLE_PKG
\endverbatim

will define the preprocessor variable <tt>USE_EXAMPLE_PKG</tt>
to the C++ compiler flags.  The codes for multiple components should
be separated by spaces, as in

\verbatim
USE_PACKAGE PKG1 PKG2
\endverbatim

Some of the supplied example applications have MPI support, which is 
activated with

\verbatim
USE_PACKAGE MPI
\endverbatim

resulting in the preprocessor definition <tt>USE_MPI</tt>.


\subsubsection subsubsec_build_imagemagick ImageMagick Functionality

The %QuickFlash library includes optional support for image output 
via the <a href="http://www.imagemagick.org">ImageMagick</a> C++ bindings (Magick++).  To activate image
support, add the code <tt>MAGICK</tt> to the <tt>USE_PACKAGE</tt>
setting in <tt>machines.cfg</tt>.

When using ImageMagick, %QuickFlash will expect the following 
definitions to be set in <tt>machines.cfg</tt>:

<ul>
  <li><tt>CCFLAGS_MAGICK</tt> &mdash; Compiler flags for finding the ImageMagick header files</li>
  <li><tt>LDFLAGS_MAGICK</tt> &mdash; Linker flags for binding to the ImageMagick library</li>
</ul>

\subsection subsec_comp_build_py Using the Build System

The <tt>build.py</tt> script performs several tasks, which are controlled by the command line flags that are provided.  The most important are the following:

<ul>
  <li><tt>--help</tt> &mdash; The complete list of options</li>
  <li><tt>--machine=</tt><em>context_name</em> &mdash; Specifies the name of the context in <tt>machines.cfg</tt> to be used; if this option is omitted, a new <tt>settings.make</tt> will not be generated, and the existing one will be used for building the library</li>
  <li><tt>--machine-config=</tt><em>config_file</em> &mdash; Use a configuration file other than <tt>machines.cfg</tt></li>
  <li><tt>--machine-help</tt> &mdash; Prints the names and descriptions of contexts defined in <tt>machines.cfg</tt></li>
  <li><tt>--screen-only</tt> &mdash; Sends the contents of the generated <tt>settings.make</tt> file to standard output; does <em>not</em> build the library or change any files</li>
  <li><tt>--enable-python</tt> &mdash; Builds both the %QuickFlash library and the Python bindings</li>
  <li><tt>--python-only</tt> &mdash; Build (or clean) only the Python bindings; the C++ library must have already been built</li>
  <li><tt>--clean</tt> &mdash; Cleans the build directories; must be used with <tt>--enable-python</tt> to clean the <tt>swig</tt> directory</li>
  <li><tt>--clean-first</tt> &mdash; Perform a clean operation before any specified or implide build operations</li>
  <li><tt>--make-cmd=</tt><em>make_command</em> &mdash; Specifies the make command to be used; overrides the <tt>MAKE</tt> setting in <tt>machines.cfg</tt></li>
  <li><tt>--makefile-only</tt> &mdash; Generates <tt>settings.make</tt> but does not build the library</li>
</ul>

Running <tt>build.py</tt> without arguments will result in an attempt
to build the C++ library using an existing <tt>settings.make</tt> file. 
This file can be one created by previous calls to <tt>build.py</tt>,
or one created by the user.  (User-created <tt>settings.make</tt>
files are generally only used for debugging purposes, as
subsequent <tt>build.py</tt> calls using the <tt>--machines=...</tt>
argument will overwrite any existing <tt>settings.make</tt> file.)


*/
