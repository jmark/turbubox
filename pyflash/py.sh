#!/usr/bin/env bash

export SCRIPT_DIR="$(dirname "${BASH_SOURCE[0]}")"

if [[ "$(hostname)" =~ 'cheops' ]]
then
    module load hdf5/1.8.13     2> /dev/null
    module load openmpi         2> /dev/null
    module load python/3.4.3    2> /dev/null
fi

export PYTHONPATH="$PYTHONPATH:$SCRIPT_DIR"
exec python3 -u $*
