#!/bin/env gnuplot

set terminal svg enhanced dashed
set encoding utf8

set size ratio 0.6

mach		= `echo $GPI_mach`
crosstime	= `echo $GPI_crosstime`
datafiles	= "`echo $GPI_datafiles`"
legends		= "`echo $GPI_legends`"

set logscale xy
set format y '%0.0e'

scale_x(x) = x
scale_y(x) = x

set grid
set xrange [1:130]
set xtics 1,2,128

set xlabel 'grid number {/:Italic k}   →'
set ylabel 'power transfer |{/:Italic P}|   →'

plot for [i=1:words(datafiles)] word(datafiles,i) \
		u (scale_x($1)):(scale_y($2)) w l lw 2 t word(legends,i)
