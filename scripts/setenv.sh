#!/bin/sh

export SCRIPT_NAME="$(basename  "$0")"
export SCRIPT_DIR="$(readlink -m "$0")"
export PROJECT_DIR="${SCRIPT_DIR}/.."

if expr "$(hostname)" : '^cheops' > /dev/null
then
    module load hdf5/1.8.13     2> /dev/null
    module load openmpi         2> /dev/null
    module load python/3.4.3    2> /dev/null
fi

export PYTHONPATH="$PYTHONPATH:$PROJECT_DIR/lib:$PROJECT_DIR/bin"
