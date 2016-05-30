#!/usr/bin/env bash

# SPATH="${BASH_SOURCE[0]}"
# PPATH="$(dirname "${BASH_SOURCE[0]}")"

PPATH="/srv/projects/astro/toolbox/turbubox"

#PCMD="~/frameworks/parallel/bin/parallel"
PCMD="parallel -k"

DIR="$1"

SRC_DIR="$DIR/hdf5"

LOG="$DIR/log/parallel-hdf5-$(date +%Y-%m-%d)"

#mkdir -p "$OUT"
#mkdir -p "$(dirname "$LOG")"

CMD=$(cat <<end
echo -n '{#}';
echo -ne '\t';
$PPATH/bin/time-evolution {}
end
)

find $SRC_DIR -type f | $PCMD --joblog $LOG $DRY $CMD
