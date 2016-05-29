include env.mk/$(shell hostname).mk

# disable multithreading
CCFLAGS_FFTW += -DFFTWPP_SINGLE_THREAD

# disable optimzations for faster compilation
#CCFLAGS += -O0

SRC_DIR = src
BIN_DIR = bin

qfl2hdf5:                   $(BIN_DIR)/qfl2hdf5
slice:						$(BIN_DIR)/slice
spectrum:                   $(BIN_DIR)/spectrum
powerspectrum:              $(BIN_DIR)/powerspectrum
spherical-shell-spectrum:   $(BIN_DIR)/spherical-shell-spectrum
root-mean-square:           $(BIN_DIR)/root-mean-square
box-average:                $(BIN_DIR)/box-average

clean:
	@rm -vf bin/* 

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

$(BIN_DIR)/slice: $(SRC_DIR)/slice.cpp
	$(CXX) -o $@ $< \
        $(CCFLAGS) \
        $(LDFLAGS) \
        -I${FFTWPP} \
        -I${HIGHFIVE} \
        ${CCFLAGS_HDF5} \
        ${LDFLAGS_HDF5}

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

$(BIN_DIR)/powerspectrum: $(SRC_DIR)/powerspectrum.cpp
	$(CXX) -o $@ $< \
        $(CCFLAGS) \
        $(LDFLAGS) \
        ${FFTWPP}/fftw++.cc \
        -I${FFTWPP} \
        -I${HIGHFIVE} \
        ${CCFLAGS_HDF5} \
        ${LDFLAGS_HDF5} \
		$(LDFLAGS_FFTW)

$(BIN_DIR)/spherical-shell-spectrum: $(SRC_DIR)/spherical-shell-spectrum.cpp
	$(CXX) -o $@ $< \
        $(CCFLAGS) \
        $(LDFLAGS) \
        ${FFTWPP}/fftw++.cc \
        -I${FFTWPP} \
        -I${HIGHFIVE} \
        ${CCFLAGS_HDF5} \
        ${LDFLAGS_HDF5} \
		$(LDFLAGS_FFTW)

$(BIN_DIR)/box-average: $(SRC_DIR)/box-average.cpp
	$(CXX) -o $@ $< \
        $(CCFLAGS) \
        $(LDFLAGS) \
        -I${FFTWPP} \
        -I${HIGHFIVE} \
        ${CCFLAGS_HDF5} \
        ${LDFLAGS_HDF5}

$(BIN_DIR)/root-mean-square: $(SRC_DIR)/root-mean-square.cpp
	$(CXX) -o $@ $< \
        $(CCFLAGS) \
        $(LDFLAGS) \
        -I${FFTWPP} \
        -I${HIGHFIVE} \
        ${CCFLAGS_HDF5} \
        ${LDFLAGS_HDF5}
