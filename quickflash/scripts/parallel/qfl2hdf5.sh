#!/usr/bin/env bash

export PRJ_DIR="$(dirname "${BASH_SOURCE[0]}")/../.."
export SRC_DIR="${1:?No source directory given!}"
export OUT_DIR="${2:?No output directory given!}"
export LOG_FIL="${3:?No logging directory given!}"

if [ -v DEBUG ]
then
    export DBG_ECHO="echo"
fi

if [ -v DRYRUN ]
then
    export DRYRUN="--dry-run"
fi

export DATASETS="dens velx vely velz accx accy accz magx magy magz"

process()
{
    local SRC_FILE="$1"
    local OUT_FILE="$OUT_DIR/$(basename $SRC_FILE)"
    
    $DBG_ECHO $PRJ_DIR/bin/qfl2hdf5 $SRC_FILE $OUT_FILE $DATASETS \
    && echo "$SRC_FILE finnished!"
}

readonly -f process
export -f process

find $SRC_DIR -type f \
    | grep -P '/\d{4}$' \
    | sort \
    | parallel --resume --joblog $LOG_FIL $DRYRUN process
