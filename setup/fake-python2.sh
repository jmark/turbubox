#!/usr/bin/env bash

if echo "$1" | grep -qE 'help|usage'
then
	echo "usage: $(basename $0) (no parameters)"
	exit 1
fi

set -eu

if expr "$(hostname)" : '^jmark' > /dev/null
then
	export ENV_CFLAGS_HDF5=-I/usr/include/hdf5_18
	export ENV_LIB_HDF5=-L/usr/lib/hdf5_18

	#export ENV_CFLAGS_HDF5=-I$HOME/builds/bin/hdf5-1.6.9/include
	#export ENV_LIB_HDF5=-L$HOME/builds/bin/hdf5-1.6.9/lib

    TMPDIR=$(mktemp -d)
    trap "rm -r $TMPDIR" EXIT

    # mask possible python3 default installation
    ln -sf $(which python2) $TMPDIR/python 
    export PATH="/$TMPDIR:$PATH"
fi

eval "$@"
