#!/bin/sh

set -e

module purge
module load gnu/5.1.0
module load openmpi/2.1.1

OMPI="/opt/rrzk/lib/openmpi/2.1.1/gcc/bin"

export CC="${OMPI}/mpicc"
export CXX="${OMPI}/mpicc"
export FC="${OMPI}/mpifort"
export F9X="${OMPI}/mpif90"
export MPI_C_COMPILER="${OMPI}/mpicc"
export MPI_CXX_COMPILER="${OMPI}/mpicc"
export MPI_Fortran_COMPILER="${OMPI}/mpifort"

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/home/jmarker2/builds/bin/hdf5-1.10.1/lib"

exec "$@"
