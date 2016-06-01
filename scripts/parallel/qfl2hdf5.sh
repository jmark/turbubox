#!/usr/bin/env bash

PROJECT_DIR="/srv/projects/astro/toolbox/turbubox"

#PARALLEL_CMD="~/frameworks/parallel/bin/parallel"
PARALLEL_CMD="parallel"

PWD_DIR="${1:?No working directory given!}"
SRC_DIR="$PWD_DIR/amr" # adaptive mesh refinement (FLASH file format)
OUT_DIR="$PWD_DIR/ugm" # uniform grid mesh (my file format)
LOG_DIR="$PWD_DIR/log" # logging information
LOG_FILE="${LOG_DIR}/qfl2hdf5_$(date +%Y-%m-%d-%H-%M-%S)"

if [ -v DRYRUN ]
then
    DRYRUN="--dry-run"
fi

DATASETS="dens velx vely velz accx accy accz magx magy magz"

CMD=$(cat <<end
${PROJECT_DIR}/bin/qfl2hdf5 {} $OUT_DIR/{/} $DATASETS
&& echo '{/} finnished!'
end
)

mkdir -p "$OUT_DIR"
mkdir -p "$LOG_DIR"

find $SRC_DIR/flash_hdf5_plt_cnt_* -type f \
| ${PARALLEL_CMD} --joblog $LOG_FILE $DRYRUN $CMD
