#!/usr/bin/env bash

module load hdf5
module load intel
module load intelmpi

pushd "$1"
    make -j4
popd
