#!/usr/bin/env bash

export PRJ_DIR="$(dirname "${BASH_SOURCE[0]}")/../.."
export SRC_DIR="${1:?No source directory given!}"

if [ -v DEBUG ]
then
    export DBG_ECHO="echo"
fi

if [ -v DRYRUN ]
then
    export DRYRUN="--dry-run"
fi

process()
{
    local SRC_FILE="$1"
    $DBG_ECHO $PRJ_DIR/bin/time-evolution $SRC_FILE
}

readonly -f process
export -f process

find $SRC_DIR -type f \
    | grep -P '/\d{4}$' \
    | sort \
    | parallel -k $DRYRUN process
