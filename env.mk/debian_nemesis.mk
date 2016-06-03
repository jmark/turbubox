CXX = g++

QUICKFLASH  = ../../frameworks/QuickFlash-1.0.0
FFTWPP      = ../../frameworks/fftw++-2.02
HIGHFIVE    = ../../frameworks/HighFive/include/highfive

CCFLAGS_HDF5 = -DH5_USE_16_API -I/usr/include/hdf5/mpich -I/usr/include/mpich
CCFLAGS_QUFL = -I$(QUICKFLASH)/include
CCFLAGS_FFTW = -I$(FFTWP)
CCFLAGS_HIGHFIVE = -I$(HIGHFIVE)

CCFLAGS += -std=c++11
CCFLAGS += -Wall
CCFLAGS += -O3

LDFLAGS_FFTW = -L/usr/lib/x86_64-linux-gnu -lfftw3
LDFLAGS_HDF5 = -L/usr/lib/x86_64-linux-gnu/hdf5/mpich -lhdf5
LDFLAGS_QUFL = $(QUICKFLASH)/lib/libquickflash.so
