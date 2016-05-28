include env.mk/$(shell hostname).mk

# disable multithreading
CCFLAGS_FFTW += -DFFTWPP_SINGLE_THREAD

# disable optimzations for faster compilation
CCFLAGS += -O0

SRC_DIR = src
BIN_DIR = bin

qfl2hdf5:       $(BIN_DIR)/qfl2hdf5
z-projection:   $(BIN_DIR)/z-projection
powerspectrum:  $(BIN_DIR)/powerspectrum
qfltest:        $(BIN_DIR)/qfltest
spectrum:       $(BIN_DIR)/spectrum

$(BIN_DIR)/qfl2hdf5: $(SRC_DIR)/qfl2hdf5.cpp
	$(CXX) -o $@ $< \
        ${CCFLAGS} \
        ${LDFLAGS} \
        -I${FFTWPP} \
        -I${HIGHFIVE} \
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
        -I${FFTWPP} \
        -I${HIGHFIVE} \
        $(CCFLAGS) \
        $(LDFLAGS)

$(BIN_DIR)/spectrum: $(SRC_DIR)/spectrum.cpp
	$(CXX) -o $@ $< \
        $(CCFLAGS) \
        $(LDFLAGS) \
        ${FFTWPP}/fftw++.cc \
        -I${FFTWPP} \
        -I${HIGHFIVE} \
        ${CCFLAGS_HDF5} \
        ${LDFLAGS_HDF5} \
		$(LDFLAGS_FFTW)

$(BIN_DIR)/qfltest: $(SRC_DIR)/qfltest.cpp
	$(CXX) -o $@ $< \
        $(CCFLAGS) \
        ${CCFLAGS_HDF5} \
        ${LDFLAGS_HDF5} \
        ${CCFLAGS_QUFL} \
        ${LDFLAGS_QUFL}
