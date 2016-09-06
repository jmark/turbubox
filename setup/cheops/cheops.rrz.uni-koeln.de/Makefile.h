FCOMP   = mpiifort
CCOMP   = mpiicc
CPPCOMP = mpiicpc
LINK    = mpiifort

PP      = -D

FFLAGS_OPT   = -c -r8 -i4 -O3 -real_size 64
FFLAGS_DEBUG = -c -g -r8 -i4 -check bounds \
					-check format -check output_conversion \
					-warn all -check uninit -traceback \
					-fp-stack-check -real_size 64
FFLAGS_TEST  = -c -r8 -i4 -O2 -real_size 64

CFLAGS_OPT   = -c -O3 -D_LARGEFILE64_SOURCE
CFLAGS_DEBUG = -c -g -debug extended -D_LARGEFILE64_SOURCE
CFLAGS_TEST  = -c -O2 -D_LARGEFILE64_SOURCE

CFLAGS_HDF5  = -DH5_USE_16_API
CFLAGS_NCMPI =
CFLAGS_MPI   =

LFLAGS_OPT   = -r8 -i4 -Ur -o
LFLAGS_DEBUG = -r8 -i4 -g -o
LFLAGS_TEST  = -r8 -i4 -o

LIB_OPT   = 
LIB_DEBUG = 
LIB_TEST  =

LIB_HDF4  = 
LIB_HDF5  = -lhdf5 -lz -lstdc++

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
