# Makefile for setup [ nemesis ]
# My supi-dupi Arch Linux setup

SWIG = swig
CCFLAGS_GSL =  
CCFLAGS_MAGICK = -I/usr/include/ImageMagick-6
ARFLAGS = rsc
USE_PACKAGE = MAGICK
CCFLAGS_NOSHARED = -mdynamic-no-pic
LD = g++
CCFLAGS_HDF5 = -DH5_USE_16_API
LDFLAGS = -W -Wall -g -O3 -ftree-vectorize -falign-loops=16
SO_EXT = so
LD_PYTHON = g++
LDFLAGS_NOSHARED = -mdynamic-no-pic
CCFLAGS_SHARED_PYTHON = -fPIC
LDFLAGS_SHARED = -shared
LDFLAGS_PYTHON =  
PYTHON_INSTALL_DIR = ./python
RANLIB = ranlib
AR = ar
CXX_PYTHON = g++
CCFLAGS_PYTHON = -I/usr/include/python2.7
CXX = g++
LDFLAGS_GSL = -lgsl -lgslcblas
MAKE = make -j4
LDFLAGS_HDF5 = -lhdf5
LDFLAGS_SHARED_PYTHON = -shared
CCFLAGS_SHARED = -fPIC
CCFLAGS = -W -Wall -g -O3 -ftree-vectorize -falign-loops=16 -DUSE_MAGICK
LIB_INSTALL_DIR = .
LDFLAGS_MAGICK = -lMagick++-6.Q16HDRI
