# Universal Makefile for FLASH Code 4.3

CCOMP   = mpiicc
CPPCOMP = mpiicpc

FCOMP   = mpiifort
LINK    = mpiifort

PP      = -D

FFLAGS_OPT   = -c -r8 -i4 -O3 -real_size 64
FFLAGS_TEST  = -c -r8 -i4 -O2 -real_size 64
FFLAGS_DEBUG = -c -r8 -i4 -O1 -real_size 64 \
				-g -check bounds -check format \
				-warn all -check uninit -traceback \
				-fp-stack-check -check output_conversion

CFLAGS_OPT   = -c -O3 -D_LARGEFILE64_SOURCE
CFLAGS_TEST  = -c -O2 -D_LARGEFILE64_SOURCE
CFLAGS_DEBUG = -c -O1 -D_LARGEFILE64_SOURCE -g -debug extended 

CFLAGS_HDF5  = $(ENV_CFLAGS_HDF5) -DH5_USE_16_API
CFLAGS_NCMPI =
CFLAGS_MPI   =

LFLAGS_OPT   = -r8 -i4 -o
LFLAGS_TEST  = -r8 -i4 -o
LFLAGS_DEBUG = -r8 -i4 -g -o

LIB_OPT   = 
LIB_DEBUG = 
LIB_TEST  =

LIB_HDF4  = 
LIB_HDF5  = $(ENV_LIB_HDF5) -lhdf5 -lz -lstdc++

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
