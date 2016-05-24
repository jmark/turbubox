# Makefile for setup [ debian-FLASH4 ]
# systemd-nspawn debian

SWIG = swig
CCFLAGS_GSL = -I /usr/include
CCFLAGS_MAGICK =  
ARFLAGS = rsc
USE_PACKAGE =  
CCFLAGS_NOSHARED = -mdynamic-no-pic
LD = g++
CCFLAGS_HDF5 = -DH5_USE_16_API -I/usr/include/hdf5/mpich -I/usr/include/mpich
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
CCFLAGS_PYTHON = -I /usr/include/python
CXX = g++
LDFLAGS_GSL = -L/usr/lib -lgsl -lgslcblas
MAKE = make -j16
LDFLAGS_HDF5 = -L/usr/lib/x86_64-linux-gnu/hdf5/mpich -lhdf5
LDFLAGS_SHARED_PYTHON = -shared
CCFLAGS_SHARED = -fPIC
CCFLAGS = -W -Wall -g -O3 -ftree-vectorize -falign-loops=16
LIB_INSTALL_DIR = .
LDFLAGS_MAGICK =  
