#!/usr/bin/env bash

if echo "$1" | grep -qE 'help|usage'
then
	echo "usage: $(basename $0) (no parameters)"
	exit 1
fi

set -eu

if expr "$(hostname)" : '^cheops' > /dev/null
then
    module purge
    module load gnu/4.8.2
    module load openmpi/1.8.6
    module load hdf5
    module load cmake

elif expr "$(hostname)" : '^jmark' > /dev/null
then
	export ENV_CFLAGS_HDF5=-I/usr/include/hdf5_18
	export ENV_LIB_HDF5=-L/usr/lib/hdf5_18

    TMPDIR=$(mktemp -d)
    trap "rm -r $TMPDIR" EXIT

    # mask possible python3 default installation
    ln -sf $(which python2) $TMPDIR/python 
    export PATH="/$TMPDIR:$PATH"
fi

eval "$@"
