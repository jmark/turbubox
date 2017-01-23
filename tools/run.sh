#!/bin/sh


# load CHEOPS modules only if we're on it.
if hostname | grep -q '^cheops'
then
    # Hardcode project directory.
    export PROJECT_DIR="$HOME/turbubox/tools"
    module purge
    module load hdf5/1.8.13     2> /dev/null
    module load openmpi         2> /dev/null
    module load python/3.4.3    2> /dev/null

elif hostname | grep -q 'jmark'
then
    export PROJECT_DIR="$HOME/projects/stirturb/turbubox/tools"
fi

export PYTHONPATH="$PYTHONPATH:$PROJECT_DIR/lib:$PROJECT_DIR/bin"
eval "$@"
