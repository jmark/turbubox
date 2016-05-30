#!/usr/bin/env bash

INFILE="$1"
OUTFILE="$2"

gnuplot <<END
set terminal pngcairo
set output '$OUTFILE'

unset key

splot '$INFILE' using 4:5:7 w l
END
