#!/usr/bin/env bash

INFILE="$1"
OUTFILE="$2"

NR=$(basename "$INFILE" | perl -F_ -e 'print $F[-1]')

gnuplot <<END
set terminal pngcairo
set output '$OUTFILE'

set title 'MHD Gir. Turb. 64x64 at z = 0, n = $NR'

set xrange [0:1]
set yrange [0:1]
set zrange [0:20]

unset key

splot '$INFILE' using 4:5:7 w l
END
