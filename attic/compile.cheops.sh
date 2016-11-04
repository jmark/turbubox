#!/usr/bin/env bash

module load hdf5/1.8.13
module load openmpi
module load python/3.4.3

cython shell_avg.pyx

gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing \
    -I/opt/rrzk/software/python/Python-3.4.3/include/python3.4m \
    -o shell_avg.so shell_avg.c
