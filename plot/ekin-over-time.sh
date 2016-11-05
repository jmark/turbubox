#!/bin/bash

MACH=$1; shift
TURN=$1; shift
DATA=$1; shift
LEGE=$1; shift

gnuplot <<END
set terminal svg enhanced dashed
set encoding utf8

set size ratio 0.6

mach		= $MACH 
crosstime	= $TURN 
datafiles	= "$DATA"
legends		= "$LEGE"

scale_x(x) = x / crosstime
scale_y(x) = x

set grid

set xrange[0:5]

set xtics 0.5

set xlabel 'crossing time {/:Italic t / t_c}   →'
set ylabel 'kinetic energy {/:Italic 𝓚}   →'

plot for [i=1:words(datafiles)] word(datafiles,i) \
		u (scale_x(\$3)):(scale_y(\$4)) w l lw 2 t word(legends,i)
END
