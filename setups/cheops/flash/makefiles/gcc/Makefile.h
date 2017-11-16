# Universal Makefile for FLASH Code 4.3

CCOMP   = mpicc
CPPCOMP = mpicxx
FCOMP   = mpif90 -c
LINK    = mpif90 -o

PP      = -D

FFLAGS = -ffree-line-length-none -fdefault-real-8 -fdefault-double-8

FFLAGS_OPT   = $(FFLAGS) -O3
FFLAGS_TEST  = $(FFLAGS) -O2
FFLAGS_DEBUG = $(FFLAGS) -O1 \
				-g -check bounds -check format \
				-warn all -check uninit -traceback \
				-fp-stack-check -check output_conversion

CFLAGS_OPT   = -c -O3 -D_LARGEFILE64_SOURCE
CFLAGS_TEST  = -c -O2 -D_LARGEFILE64_SOURCE
CFLAGS_DEBUG = -c -O1 -D_LARGEFILE64_SOURCE -g -debug extended 

CFLAGS_HDF5  = $(ENV_CFLAGS_HDF5) -DH5_USE_16_API
CFLAGS_NCMPI =
CFLAGS_MPI   =

LFLAGS_OPT   =
LFLAGS_TEST  =
LFLAGS_DEBUG = -g

LIB_OPT   = 
LIB_DEBUG = 
LIB_TEST  =

LIB_HDF4  = 
LIB_HDF5  = $(ENV_LIB_HDF5) -L/opt/hdf5-1.8.19-openmpi/lib -lhdf5 -lz -lstdc++

LIB_PAPI  =
LIB_MATH  = 

LIB_MPI   = 
LIB_NCMPI =
LIB_MPE   =

MACHOBJ =

MV = mv -f
AR = ar -r
RM = rm -f
CD = cd
RL = ranlib
ECHO = echo
