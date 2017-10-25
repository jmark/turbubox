#!/bin/env gnuplot

set terminal svg enhanced dashed
set encoding utf8

set size ratio 0.6

mach		= `echo $GPI_mach`
crosstime	= `echo $GPI_crosstime`
datafiles	= "`echo $GPI_datafiles`"
legends		= "`echo $GPI_legends`"

scale_x(x) = x / crosstime
scale_y(x) = x

set grid

set xrange[0:5]

set xtics 0.5

set xlabel 'crossing time {/:Italic t / t_c}   →'
set ylabel 'kinetic energy {/:Italic 𝓚}   →'

plot for [i=1:words(datafiles)] word(datafiles,i) \
		u (scale_x($3)):(scale_y($4)) w l lw 2 t word(legends,i)
