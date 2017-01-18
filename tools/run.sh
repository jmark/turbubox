#!/bin/sh


# load CHEOPS modules only if we're on it.
<<<<<<< HEAD
if hostname | grep -q '^cheops'
=======
if hostname | grep -qE '^cheops'
>>>>>>> 41bfc6d777e409b93deb5fd8727c227541b17906
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
