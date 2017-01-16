#!/bin/sh

set -eu

if hostname | grep -qi 'cheops'
then
    module purge
    module load hdf5
    module load intel
    module load intelmpi

    eval srun -n "$NTASK" "$@"

elif hostname | grep -qi 'jmark'
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

    eval mpirun -n "$NTASK" "$@"
fi
