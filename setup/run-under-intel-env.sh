#!/usr/bin/env bash

if echo "$1" | grep -qE 'help|usage'
then
	echo "usage: $(basename $0) (no parameters)"
	exit 1
fi

#set -eu

if hostname | grep -q '^cheops'
then
    module purge

    export CC=icc
    export CXX=icc
    export FC=ifort
    export F9X=ifort
    export MPI_C_COMPILER=mpicc
    export MPI_CXX_COMPILER=mpicxx
    export MPI_Fortran_COMPILER=mpiifort
    export I_MPI_CC=icc
    export I_MPI_CXX=icpc
    export I_MPI_F77=ifort
    export I_MPI_F90=ifort

    module load hdf5
    module load intel
    module load intelmpi
    module load mkl
    module load cmake

# SuperMUC
elif hostname | grep -q '^login'
then
    echo 'Detected Supermuc!'

    source /etc/profile
    source /etc/profile.d/modules.sh

    export CC=icc
    export CXX=icc
    export FC=ifort
    export F9X=ifort
    export MPI_C_COMPILER=mpicc
    export MPI_CXX_COMPILER=mpicxx
    export MPI_Fortran_COMPILER=mpiifort
    export I_MPI_CC=icc
    export I_MPI_CXX=icpc
    export I_MPI_F77=ifort
    export I_MPI_F90=ifort

    #module purge
    module load admin
    module load tempdir
    module load lrztools
    module load poe
    module load intel
    module load mkl
    module load mpi.ibm
    module load hdf5/mpi
    module load cmake

elif hostname | grep -q 'jmark'
then
    INTEL_ROOT="/opt/intel/compilers_and_libraries_2017.0.098"
    set +u
        #source /etc/profile.d/intel*
        source $INTEL_ROOT/linux/bin/compilervars.sh -arch intel64 -platform linux
    set -u

    export CC=icc
    export CXX=icc
    export FC=ifort
    export F9X=ifort
    export MPI_C_COMPILER=mpicc
    export MPI_CXX_COMPILER=mpicxx
    export MPI_Fortran_COMPILER=mpiifort
    export I_MPI_CC=icc
    export I_MPI_CXX=icpc
    export I_MPI_F77=ifort
    export I_MPI_F90=ifort

    # export PCC=mpiicc
    # export PFC=mpiifort

    # export MPI_C_COMPILER=mpiicc
    # export MPI_CXX_COMPILER=mpiicxx
    # export MPI_Fortran_COMPILER=mpiifort

    # export I_MPI_CC=mpiicc
    # export I_MPI_CXX=mpiicpc
    # export I_MPI_F77=mpiifort
    # export I_MPI_F90=mpiifort

	#export ENV_CFLAGS_HDF5=-I/usr/include/hdf5_18
	#export ENV_LIB_HDF5=-L/usr/lib/hdf5_18

	# export ENV_CFLAGS_HDF5=-I$HOME/builds/bin/hdf5-1.6.9/include
	# export ENV_LIB_HDF5=-L$HOME/builds/bin/hdf5-1.6.9/lib

    # TMPDIR=$(mktemp -d)
    # trap "rm -r $TMPDIR" EXIT

    # # mask possible python3 default installation
    # ln -sf $(which python2) $TMPDIR/python 
    # export PATH="/$TMPDIR:$PATH"
fi

eval "$@"
