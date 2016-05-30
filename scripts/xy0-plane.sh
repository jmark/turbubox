#!/usr/bin/env bash

#DIR='/mnt/seagate/FLASH/stirring_test/2016-05-13_64'
#DIR='/mnt/seagate/FLASH/bouchut3/stirring_test/2016-05-14_64'
DIR='/mnt/seagate/FLASH/bouchut5/stirring_test/2016-05-15_64'

SNAPSHOTS="$DIR/out/*plt_cnt_*"
CSV="$DIR/csv"
LOG="$DIR/log/parallel"

find $SNAPSHOTS | parallel \
    --joblog $LOG \
    --progress \
    --eta \
    "./bin/z-projection {} > $CSV/{/}"
