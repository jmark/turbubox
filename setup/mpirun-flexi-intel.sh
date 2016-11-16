#!/bin/sh

FLEXI_EXE="${1:?No FLASH executable given!}" && shift
FLEXI_INI="${1:?No parameter file given!}" && shift

if expr "$(hostname)" : '^cheops' > /dev/null
then
    module purge
    module load hdf5
    module load intel
    module load intelmpi

    srun -n "$NTASK" "$FLASH_EXE" -par_file "$FLASH_PAR" "$@"

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

    source /etc/profile.d/intel*

    mpirun -n "$NTASK" "$FLEXI_EXE" "$FLEXI_INI" "$@"
fi
