#!/usr/bin/env bash

SPATH="${BASH_SOURCE[0]}"
PPATH="$(dirname "${BASH_SOURCE[0]}")"

PPATH="/srv/projects/astro/toolbox/turbubox"

#PCMD="~/frameworks/parallel/bin/parallel"
PCMD="parallel"

DBS="dens"

DIR="$1"

HD5="$DIR/hdf5/*"
OUT="$DIR/csv/xy0-plane"
LOG="$DIR/log/parallel"

#DRY="--dry-run"

CMD=$(cat <<end
$PPATH/bin/slice {} $DBS @ - - 0 > $OUT/{/} && echo '{/} finnished!'
end
)

mkdir -p "$OUT"
mkdir -p "$(dirname "$LOG")"

find $HD5 -type f | parallel \
    --joblog $LOG \
    $DRY \
    $CMD
