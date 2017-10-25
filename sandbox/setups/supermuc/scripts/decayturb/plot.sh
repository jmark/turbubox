#!/bin/bash

mkdir -p log cache png

NPROCS="${NPROCS:-4}"
NFILES="$(echo "$@" | wc -w)"
CRANGE="cdens=(-1.0,1.0), cekin=(0.02,2.0), cmach=(0.2,1.2), cvort=(-12,-7)"

$HOME/turbubox/plot/dens-ekin-mach-vort.py \
    --destdir png/ --cachedir cache/ \
    --parallel $NPROCS --ntasks $NFILES \
    --crange "$CRANGE" "$@"
