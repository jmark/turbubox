#!/usr/bin/env bash

SPATH="${BASH_SOURCE[0]}"
PPATH="$(dirname "${BASH_SOURCE[0]}")"

PPATH="/srv/projects/astro/toolbox/turbubox"

#PCMD="~/frameworks/parallel/bin/parallel"
PCMD="parallel"

DIR="$1"

IN_DIR="$DIR/csv/xy0-plane"
OUT_DIR="$DIR/png/xy0-plane"
LOG="$DIR/log/parallel"

mkdir -p "$OUT_DIR"

#DRY="--dry-run"

CMD=$(cat <<end
$PPATH/scripts/plots/surface.sh {} $OUT_DIR/{/}.png && echo '{/} finnished!'
end
)

find $IN_DIR -type f | parallel \
    --joblog $LOG \
    $DRY \
    $CMD
