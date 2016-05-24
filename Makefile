CXX = g++

QUICKFLASH  = $(HOME)/frameworks/QuickFlash-1.0.0
FFTWpp      = ${HOME}/frameworks/fftw++-2.02
HIGHFIVE    = ${HOME}/frameworks/HighFive/include/highfive

CCFLAGS_HDF5 = -DH5_USE_16_API -I/usr/include/hdf5/mpich -I/usr/include/mpich
LDFLAGS_HDF5 = -L/usr/lib/x86_64-linux-gnu/hdf5/mpich -lhdf5

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

SRC_DIR = src
BIN_DIR = bin

z-projection: $(BIN_DIR)/$@.cpp
$(BIND_DIR)/z-projection.cpp : $(SRC_DIR)/z-projection.cpp
	$(CXX) -o bin/$@ $< $(CCFLAGS) $(LDFLAGS)

powerspectrum: $(BIN_DIR)/$@.cpp
$(BIN_DIR)powerspectrum: $(SRC_DIR)/powerspectrum.cpp
	$(CXX) -o bin/$@ $< ${FFTWpp}/fftw++.cc $(CCFLAGS) $(LDFLAGS)

qfl2hdf5: $(BIN_DIR)/$@.cpp
$(BIN_DIR)/qfl2hdf5: $(SRC_DIR)/qfl2hdf5.cpp
	$(CXX) -o bin/$@ $< -I${HIGHFIVE} $(CCFLAGS) ${LDFLAGS_HDF5} ${LDFLAGS}
