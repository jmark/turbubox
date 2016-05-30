#!/usr/bin/env bash

DIR="$1"
#DIR='/mnt/seagate/FLASH/stirring_test/2016-05-13_64'
#DIR='/mnt/seagate/FLASH/bouchut3/stirring_test/2016-05-14_64'
#DIR='/mnt/seagate/FLASH/bouchut5/stirring_test/2016-05-15_64'

SPATH="${BASH_SOURCE[0]}"
PPATH=$(dirname "${BASH_SOURCE[0]}")

FLS="$DIR/csv/*plt_cnt_*"
LOG="$DIR/log/parallel-$(date +%Y-%m-%d-%H-%M-%S)"
OUT="$DIR/png"

find $FLS | parallel \
    --joblog $LOG \
    --progress \
    --eta \
    "$PPATH/pbox-plot.sh {} > $OUT/{/}.png"
