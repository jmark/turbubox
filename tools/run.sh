#!/bin/sh

LIB_PATH="${1:?No path to the modules given!}" && shift

# load CHEOPS modules only if we're on it.
if hostname | grep -q '^cheops'
then
    module purge
    module load hdf5/1.8.13     2> /dev/null
    module load openmpi         2> /dev/null
    module load python/3.4.3    2> /dev/null
fi

export PYTHONPATH="$PYTHONPATH:$LIB_PATH"
eval "$@"
