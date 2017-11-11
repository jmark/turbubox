#!/bin/sh -l

set -e

module purge

source /etc/profile
source /etc/profile.d/modules.sh

module unload mpi.intel
module load gcc/7
module load mpi.ompi/2.0/gcc
#module load hdf5/mpi/1.8.17_gcc

# OMPI="/opt/rrzk/lib/openmpi/2.1.1/gcc/bin"
# 
# export CC="${OMPI}/mpicc"
# export CXX="${OMPI}/mpicc"
# export FC="${OMPI}/mpifort"
# export F9X="${OMPI}/mpif90"
# export MPI_C_COMPILER="${OMPI}/mpicc"
# export MPI_CXX_COMPILER="${OMPI}/mpicc"
# export MPI_Fortran_COMPILER="${OMPI}/mpifort"
# 
# export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/home/jmarker2/builds/bin/hdf5-1.10.1/lib"

exec "$@"
