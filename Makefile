CXX = g++

QUICKFLASH = $(HOME)/frameworks/QuickFlash-1.0.0

FFTWpp = ${HOME}/frameworks/fftw++-2.02

# CCFLAGS_HDF5 = -DH5_USE_16_API -I/usr/include/hdf5/mpich -I/usr/include/mpich
# LDFLAGS_HDF5 = -L/usr/lib/x86_64-linux-gnu/hdf5/mpich -lhdf5

# includes
CCFLAGS += -I$(QUICKFLASH)/include
CCFLAGS += -I/usr/include/hdf5/mpich
CCFLAGS += -I/usr/include/mpich
CCFLAGS += -DH5_USE_16_API
CCFLAGS += -std=c++0x
CCFLAGS += -I${FFTWpp}

# debug flags
CCFLAGS += -W
CCFLAGS += -Wall
CCFLAGS += -Wno-unused-parameter
#CCFLAGS += -g

# opt flags
CCFLAGS += -O3
# CCFLAGS += -ftree-vectorize
# CCFLAGS += -falign-loops=16

LDFLAGS += $(QUICKFLASH)/lib/libquickflash.so
LDFLAGS += -lfftw3

all : z-projection powerspectrum

z-projection: src/z-projection.cpp
	$(CXX) -o bin/z-projection src/z-projection.cpp $(CCFLAGS) $(LDFLAGS)

powerspectrum: src/powerspectrum.cpp
	$(CXX) -o bin/powerspectrum ${FFTWpp}/fftw++.cc src/powerspectrum.cpp $(CCFLAGS) $(LDFLAGS)

clean : clean-obj
	rm -f $(BIN_FILES)

clean-obj :
	rm -f *.o
