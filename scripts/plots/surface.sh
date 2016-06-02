#!/usr/bin/env bash

SRC_FILE="${1:?No source file given!}"
OUT_FILE="${2:?No output file given!}"

#NR=$(basename "$INFILE" | perl -F_ -e 'print $F[-1]')
NR=$(basename "$SRC_FILE")

gnuplot <<END
set terminal pngcairo
set output '$OUT_FILE'

set title 'MHD Gir. Turb. 64x64 at z = 0, n = $NR'

set xrange [0:1]
set yrange [0:1]
set zrange [0:20]

unset key

splot '$SRC_FILE' using 4:5:7 w l
END
