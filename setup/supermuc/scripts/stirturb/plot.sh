#!/bin/bash

mkdir -p log cache png

NPROCS="${NPROCS:-4}"
FILES="$@"
NFILES="$(echo "$FILES" | wc -w)"
CRANGE="cdens=(-1.0,1.0), cekin=(0.02,1.0), cmach=(0.2,1.2), cvort=(-10,-7)"

$HOME/turbubox/plot/dens-ekin-mach-vort.py \
    --destdir png/ \
    --parallel $NPROCS --ntasks $NFILES --title 'Midpoint Hybrid' \
    --crange "$CRANGE" $FILES
