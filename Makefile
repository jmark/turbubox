CXX = g++

QUICKFLASH  = $(HOME)/frameworks/QuickFlash-1.0.0
FFTWpp      = ${HOME}/frameworks/fftw++-2.02
HIGHFIVE    = ${HOME}/frameworks/HighFive/include/highfive

CCFLAGS_HDF5 = -DH5_USE_16_API -I/usr/include/hdf5/mpich -I/usr/include/mpich
CCFLAGS_QUFL = -I$(QUICKFLASH)/include

CCFLAGS += -std=c++11
CCFLAGS += -Wall
CCFLAGS += -O3

#LDFLAGS_FFTW += -lfftw
LDFLAGS_HDF5 = -L/usr/lib/x86_64-linux-gnu/hdf5/mpich -lhdf5
LDFLAGS_QUFL = $(QUICKFLASH)/lib/libquickflash.so

z-projection: src/z-projection.cpp
	$(CXX) -o bin/$@ $< \
        $(CCFLAGS) \
        $(LDFLAGS)

powerspectrum: src/powerspectrum.cpp
	$(CXX) -o bin/$@ $< \
        ${FFTWpp}/fftw++.cc \
        $(CCFLAGS) \
        $(LDFLAGS)

qfltest: src/qfltest.cpp
	$(CXX) -o bin/$@ $< \
        $(CCFLAGS) \
        ${CCFLAGS_HDF5} \
        ${LDFLAGS_HDF5} \
        ${CCFLAGS_QUFL} \
        ${LDFLAGS_QUFL}

qfl2hdf5: src/qfl2hdf5.cpp
	$(CXX) -o bin/$@ $< \
        -I${FFTWpp} \
        -I${HIGHFIVE} \
        ${CCFLAGS} \
        ${CCFLAGS_HDF5} \
        ${LDFLAGS_HDF5} \
        ${CCFLAGS_QUFL} \
        ${LDFLAGS_QUFL}
