include env.mk/$(shell hostname).mk

# disable multithreading
CCFLAGS_FFTW += -DFFTWPP_SINGLE_THREAD

# disable optimzations for faster compilation
#CCFLAGS += -O0

SRC_DIR = src
BIN_DIR = bin

all: \
	qfl2hdf5 \
	slice \
	spectrum \
	powerspectrum \
	spherical-shell-spectrum \
	root-mean-square \
	box-average

qfl2hdf5:                   bin $(BIN_DIR)/qfl2hdf5
slice:                      bin $(BIN_DIR)/slice
spectrum:                   bin $(BIN_DIR)/spectrum
powerspectrum:              bin $(BIN_DIR)/powerspectrum
spherical-shell-spectrum:   bin $(BIN_DIR)/spherical-shell-spectrum
root-mean-square:           bin $(BIN_DIR)/root-mean-square
box-average:                bin $(BIN_DIR)/box-average
time-evolution:             bin $(BIN_DIR)/time-evolution

clean:
	@rm -vf bin/* 

bin:
	mkdir -p bin

$(BIN_DIR)/qfl2hdf5: $(SRC_DIR)/qfl2hdf5.cpp $(SRC_DIR)/ulz.hpp
	$(CXX) -o $@ $< \
        ${CCFLAGS} \
        ${LDFLAGS} \
        $(CCFLAGS_FFTW) \
        $(CCFLAGS_HIGHFIVE) \
        ${CCFLAGS_HDF5} \
        ${LDFLAGS_HDF5} \
        ${CCFLAGS_QUFL} \
        ${LDFLAGS_QUFL}

$(BIN_DIR)/slice: $(SRC_DIR)/slice.cpp $(SRC_DIR)/ulz.hpp
	$(CXX) -o $@ $< \
        $(CCFLAGS) \
        $(LDFLAGS) \
        $(CCFLAGS_FFTW) \
        $(CCFLAGS_HIGHFIVE) \
        ${CCFLAGS_HDF5} \
        ${LDFLAGS_HDF5}

$(BIN_DIR)/spectrum: $(SRC_DIR)/spectrum.cpp $(SRC_DIR)/ulz.hpp
	$(CXX) -o $@ $< \
        $(CCFLAGS) \
        $(LDFLAGS) \
        ${FFTWPP}/fftw++.cc \
        $(CCFLAGS_FFTW) \
        $(CCFLAGS_HIGHFIVE) \
        ${CCFLAGS_HDF5} \
        ${LDFLAGS_HDF5} \
		$(LDFLAGS_FFTW)

$(BIN_DIR)/powerspectrum: $(SRC_DIR)/powerspectrum.cpp $(SRC_DIR)/ulz.hpp
	$(CXX) -o $@ $< \
        $(CCFLAGS) \
        $(LDFLAGS) \
        ${FFTWPP}/fftw++.cc \
        $(CCFLAGS_FFTW) \
        $(CCFLAGS_HIGHFIVE) \
        ${CCFLAGS_HDF5} \
        ${LDFLAGS_HDF5} \
		$(LDFLAGS_FFTW)

$(BIN_DIR)/spherical-shell-spectrum: $(SRC_DIR)/spherical-shell-spectrum.cpp $(SRC_DIR)/ulz.hpp
	$(CXX) -o $@ $< \
        $(CCFLAGS) \
        $(LDFLAGS) \
        ${FFTWPP}/fftw++.cc \
        $(CCFLAGS_FFTW) \
        $(CCFLAGS_HIGHFIVE) \
        ${CCFLAGS_HDF5} \
        ${LDFLAGS_HDF5} \
		$(LDFLAGS_FFTW)

$(BIN_DIR)/box-average: $(SRC_DIR)/box-average.cpp $(SRC_DIR)/ulz.hpp
	$(CXX) -o $@ $< \
        $(CCFLAGS) \
        $(LDFLAGS) \
        $(CCFLAGS_FFTW) \
        $(CCFLAGS_HIGHFIVE) \
        ${CCFLAGS_HDF5} \
        ${LDFLAGS_HDF5}

$(BIN_DIR)/time-evolution: $(SRC_DIR)/time-evolution.cpp $(SRC_DIR)/ulz.hpp
	$(CXX) -o $@ $< \
        $(CCFLAGS) \
        $(LDFLAGS) \
        $(CCFLAGS_FFTW) \
        $(CCFLAGS_HIGHFIVE) \
        ${CCFLAGS_HDF5} \
        ${LDFLAGS_HDF5}

$(BIN_DIR)/root-mean-square: $(SRC_DIR)/root-mean-square.cpp $(SRC_DIR)/ulz.hpp
	$(CXX) -o $@ $< \
        $(CCFLAGS) \
        $(LDFLAGS) \
        $(CCFLAGS_FFTW) \
        $(CCFLAGS_HIGHFIVE) \
        ${CCFLAGS_HDF5} \
        ${LDFLAGS_HDF5}
