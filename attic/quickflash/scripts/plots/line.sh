#!/usr/bin/env bash

INFILE="$1"
OUTFILE="$2"

NR=$(basename "$INFILE" | perl -F_ -e 'print $F[-1]')

gnuplot <<END
set terminal pngcairo
set output '$OUTFILE'

set grid

set style line 1 lw 2 pt 3 ps 0.5
set style line 2 lw 2 pt 3 ps 0.5
set style line 3 lw 2 pt 3 ps 0.5

set title "(0,0,z) line: n = $NR"

set xrange [0: 1]
set yrange [0: 5]

norm3d(x,y,z) = x*x + y*y + z*z

plot \
    '$INFILE' using 6:( 1./10.  * \$7) w l t "density" ls 1,\
    '$INFILE' using 6:( 1./20.  * norm3d(\$8,\$9,\$10)) w l t 'energy' ls 2,\
    '$INFILE' using 6:( 10**7 * norm3d(\$11,\$12,\$13)) w l t "mag. strength" ls 3
END
