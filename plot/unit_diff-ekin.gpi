#!/bin/env gnuplot

set terminal svg enhanced dashed
set encoding utf8

set size ratio 0.6

mach		= `echo $GPI_mach`
crosstime	= `echo $GPI_crosstime`
datafiles	= "`echo $GPI_datafiles`"
legends		= "`echo $GPI_legends`"

scale_x(x) = x / crosstime
scale_y(x) = log10(abs(x))

set grid

set xtics 0.2
set xrange [2:5]

set xlabel 'crossing time {/:Italic t / t_c}   →'
set ylabel 'kinetic energy transfer rate log_{10}|{/:Italic d𝓚/dt}|   →'

dfile = 'data/unit_128_b3-b5-es_mach-10_diff-ekin.dsv'

plot \
	dfile u (scale_x($1)):(scale_y($5)) w l lw 2 t 'b3',\
	dfile u (scale_x($1)):(scale_y($6)) w l lw 2 t 'b5',\
	dfile u (scale_x($1)):(scale_y($7)) w l lw 2 t 'es'
