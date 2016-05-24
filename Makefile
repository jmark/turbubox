CXX = g++

QUICKFLASH  = $(HOME)/frameworks/QuickFlash
FFTWPP      = ${HOME}/frameworks/fftw++-2.02
HIGHFIVE    = ${HOME}/frameworks/HighFive/include/highfive

CCFLAGS_HDF5 = -DH5_USE_16_API -I/usr/include/hdf5/mpich -I/usr/include/mpich
CCFLAGS_QUFL = -I$(QUICKFLASH)/include
CCFLAGS_FFTW = -I$(FFTWP)

CCFLAGS += -std=c++11
CCFLAGS += -Wall
CCFLAGS += -O3

LDFLAGS_FFTW = -lfftw
LDFLAGS_HDF5 = -L/usr/lib/x86_64-linux-gnu/hdf5/mpich -lhdf5
LDFLAGS_QUFL = $(QUICKFLASH)/bin/libquickflash.so

SRC_DIR = src
BIN_DIR = bin

qfl2hdf5:       $(BIN_DIR)/qfl2hdf5
z-projection:   $(BIN_DIR)/z-projection
powerspectrum:  $(BIN_DIR)/powerspectrum
qfltest:        $(BIN_DIR)/qfltest

$(BIN_DIR)/qfl2hdf5: $(SRC_DIR)/qfl2hdf5.cpp
	$(CXX) -o $@ $< \
        -I${FFTWPP} \
        -I${HIGHFIVE} \
        ${CCFLAGS} \
        ${CCFLAGS_HDF5} \
        ${LDFLAGS_HDF5} \
        ${CCFLAGS_QUFL} \
        ${LDFLAGS_QUFL}

$(BIN_DIR)/z-projection: $(SRC_DIR)/z-projection.cpp
	$(CXX) -o $@ $< \
        $(CCFLAGS) \
        $(LDFLAGS) \
        ${CCFLAGS_HDF5} \
        ${LDFLAGS_HDF5} \
        ${CCFLAGS_QUFL} \
        ${LDFLAGS_QUFL}

$(BIN_DIR)/powerspectrum: $(SRC_DIR)/powerspectrum.cpp
	$(CXX) -o $@ $< \
        ${FFTWPP}/fftw++.cc \
        $(CCFLAGS) \
        $(LDFLAGS)

$(BIN_DIR)/qfltest: $(SRC_DIR)/qfltest.cpp
	$(CXX) -o $@ $< \
        $(CCFLAGS) \
        ${CCFLAGS_HDF5} \
        ${LDFLAGS_HDF5} \
        ${CCFLAGS_QUFL} \
        ${LDFLAGS_QUFL}
