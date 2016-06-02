#!/usr/bin/env bash

PRJ_DIR="$(dirname "${BASH_SOURCE[0]}")/../.."

#PARALLEL_CMD="~/frameworks/parallel/bin/parallel"
PARALLEL_CMD="parallel"

PWD_DIR="${1:?No working directory given!}"
SRC_DIR="$PWD_DIR/ugm" # uniform grid mesh (my file format)
LOG_DIR="$PWD_DIR/log" # logging information
LOG_FILE="${LOG_DIR}/time-evolution_$(date +%Y-%m-%d-%H-%M-%S)"

if [ -v DRYRUN ]
then
    DRYRUN="--dry-run"
fi

mkdir -p "$LOG_DIR"

CMD=$(cat <<end
echo -n '{#}';
echo -ne '\t';
$PRJ_DIR/bin/time-evolution {}
end
)

find $SRC_DIR/* -type f \
| ${PARALLEL_CMD} -k --joblog $LOG_FILE $DRYRUN $CMD
