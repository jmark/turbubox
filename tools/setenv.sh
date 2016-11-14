#!/bin/sh

# This scripts loads necessary modules and sets environment variables in order
# to ensure proper functionality.
#
# Copying and adapting this file to your own needs is highly encouraged.
#
# Following two modes of operation are presented.
#
# 1.) interactive session
# -----------------------
# In an interactive shell session one has to 'source'
# this script like so:
#
#   $ source setenv.sh ${project_path}/tools
# 
# Then run the programs like normal.
#
# 2.) preloader script
# --------------------
# Especially interesting in batch processing w/o user intervention. See
# comments below.

# ======
# Mode 1
# ======

export PROJECT_DIR=$(readlink -m "${1:? No project directory given!}")

if expr "$(hostname)" : '^cheops' > /dev/null
then
    module purge
    module load hdf5/1.8.13
    module load openmpi
    module load python/3.4.3
fi

export PYTHONPATH="$PYTHONPATH:$PROJECT_DIR/lib:$PROJECT_DIR/bin"

# ======
# Mode 2
# ======
if false # Mode 2 does not get run per default.
then

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

# Run script file with parameters.
python3 "$@"

fi
