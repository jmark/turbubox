#!/usr/bin/env bash

SPATH="${BASH_SOURCE[0]}"
PPATH="$(dirname "${BASH_SOURCE[0]}")"

PARALLEL="~/frameworks/parallel/bin/parallel"

DIR="$1"

SPT="$DIR/out/snapshot_hdf5_plt_cnt_*"
OUT="$DIR/hdf5"
LOG="$DIR/log/parallel-hdf5-$(date +%Y-%m-%d)"

# bouchut
# DBS="dens temp pres velx vely velz accx accy accz magx magy magz magp"

# 8wave
DBS="dens velx vely velz accx accy accz magx magy magz"

rm -f "$LOG"
mkdir -p "$OUT"

find $SPT | $PARALLEL \
    --joblog $LOG \
    "$PPATH/../bin/qfl2hdf5 {} $OUT/{/} $DBS && echo '{/} finnished!'"
