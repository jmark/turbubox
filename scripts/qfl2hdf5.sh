#!/usr/bin/env bash

#DIR='/mnt/seagate/FLASH/stirring_test/2016-05-13_64'
#DIR='/mnt/seagate/FLASH/bouchut3/stirring_test/2016-05-14_64'
DIR='/mnt/seagate/FLASH/bouchut5/stirring_test/2016-05-15_64'

SNAPSHOTS="$DIR/out/snapshot_hdf5_plt_cnt_*"
OUT="$DIR/hdf5"
LOG="$DIR/log/parallel-hdf5"

DBS="dens temp pres velx vely velz accx accy accz magx magy magz magp"

rm -f "$LOG"
mkdir -p "$OUT"

find $SNAPSHOTS | parallel \
    --joblog $LOG \
    --progress \
    --eta \
    "./bin/qfl2hdf5 {} $OUT/{/} $DBS"
