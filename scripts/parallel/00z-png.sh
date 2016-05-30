#!/usr/bin/env bash

PPATH="/srv/projects/astro/toolbox/turbubox"

#PCMD="~/frameworks/parallel/bin/parallel"
PCMD="parallel"

DIR="$1"

IN_DIR="$DIR/csv/00z"
OUT_DIR="$DIR/png/00z"
LOG="$DIR/log/parallel"

mkdir -p "$OUT_DIR"

#DRY="--dry-run"

CMD=$(cat <<end
$PPATH/scripts/plots/line.sh {} $OUT_DIR/{/}.png && echo '{/} finnished!'
end
)

find $IN_DIR -type f | $PCMD --joblog $LOG $DRY $CMD
