#!/usr/bin/env bash

export PRJ_DIR="$(dirname "${BASH_SOURCE[0]}")/../.."

export SRC_DIR="${1:?No source directory given!}"
export OUT_DIR="${2:?No output directory given!}"
export AXE_DEF="${3:?No axes definition given!}"

mkdir -v -p "$OUT_DIR"

if [ -v DEBUG ]
then
    export DBG_ECHO="echo"
fi

if [ -v DRYRUN ]
then
    export DRYRUN="--dry-run"
fi

#export DATASETS="dens velx vely velz accx accy accz magx magy magz"
export DATASETS="dens ekin emag"

process()
{
    local SRC_FILE="$1"
    local OUT_FILE="$OUT_DIR/$(basename $SRC_FILE)"
    
    $DBG_ECHO $PRJ_DIR/bin/slice $SRC_FILE $DATASETS @ $AXE_DEF \
    > $OUT_DIR/$OUT_FILE \
    && echo "$SRC_FILE finnished!"
}

readonly -f process
export -f process

find $SRC_DIR -type f \
    | grep -P '/\d{4}$' \
    | sort \
    | parallel $DRYRUN process
