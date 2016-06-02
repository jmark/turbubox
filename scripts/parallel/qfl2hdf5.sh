#!/usr/bin/env bash

export PRJ_DIR="$(dirname "${BASH_SOURCE[0]}")/../.."
export PWD_DIR="${1:?No working directory given!}"
export SRC_DIR="$PWD_DIR/amr" # adaptive mesh refinement (FLASH file format)
export OUT_DIR="$PWD_DIR/ugm" # uniform grid mesh (my file format)
export LOG_DIR="$PWD_DIR/log" # logging information
export LOG_FIL="${LOG_DIR}/qfl2hdf5_$(date +%Y-%m-%d-%H-%M-%S)"

if [ -v DRYRUN ]
then
    export DRYRUN="--dry-run"
else
    mkdir -v -p "$OUT_DIR"
    mkdir -v -p "$LOG_DIR"
fi

export DATASETS="dens velx vely velz accx accy accz magx magy magz"

process()
{
    local SRC_FILE="$1"
    local OUT_FILE="$OUT_DIR/$(basename $SRC_FILE)"
    
    $PRJ_DIR/bin/qfl2hdf5 $SRC_FILE $OUT_FILE $DATASETS \
    && echo "$SRC_FILE finnished!"
}

readonly -f process
export -f process

find $SRC_DIR -type f \
    | grep -P '/\d{4}$' \
    | sort \
    | parallel --joblog $LOG_FIL $DRYRUN process
