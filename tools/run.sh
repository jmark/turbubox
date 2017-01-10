#!/bin/sh

# Hardcode project directory.
export PROJECT_DIR="$HOME/turbubox/tools"

# load CHEOPS modules only if we're on it.
if expr "$(hostname)" : '^cheops' > /dev/null
then
    module purge
    module load hdf5/1.8.13     2> /dev/null
    module load openmpi         2> /dev/null
    module load python/3.4.3    2> /dev/null
fi

export PYTHONPATH="$PYTHONPATH:$PROJECT_DIR/lib:$PROJECT_DIR/bin"

eval "$@"
