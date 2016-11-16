#!/usr/bin/env bash

set -eu

if expr "$(hostname)" : '^cheops' > /dev/null
then
    module purge
    module load hdf5
    module load intel
    module load intelmpi

    make -j "$@"

elif expr "$(hostname)" : '^jmark' > /dev/null
then
    export CC=icc
    export CXX=icc
    export FC=ifort
    export MPI_C_COMPILER=mpicc
    export MPI_CXX_COMPILER=mpicxx
    export MPI_Fortran_COMPILER=mpiifort
    export I_MPI_CC=icc
    export I_MPI_CXX=icpc
    export I_MPI_F77=ifort
    export I_MPI_F90=ifort

    set +u ; source /etc/profile.d/intel* ; set -u

	export ENV_CFLAGS_HDF5=-I/usr/include/hdf5_18
	export ENV_LIB_HDF5=-L/usr/lib/hdf5_18

    TMPDIR=$(mktemp -d)
    trap "rm -r $TMPDIR" EXIT

    # mask possible python3 default installation
    ln -sf $(which python2) $TMPDIR/python 
    export PATH="/$TMPDIR:$PATH"

    make -j "$@"
fi
