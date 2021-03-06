#!/bin/env gnuplot

set terminal svg enhanced
set encoding utf8

set size ratio 0.6

datafiles	= "`echo $GPI_datafiles`"
legends		= "`echo $GPI_legends`"

mach		= `echo $GPI_mach`
crosstime	= `echo $GPI_crosstime`
domainsize  = `echo $GPI_domainsize`
ambientdens = `echo $GPI_ambientdens`
soundspeed  = `echo $GPI_soundspeed`

totalmass   = domainsize**3 * ambientdens

scale_x(x) = x / crosstime
scale_y(x) = sqrt(2.0 * x / totalmass) / soundspeed

set grid

set xrange [0:5]
set xtics 0.5

set xlabel 'crossing time: {/:Italic t / t_c}   →'
set ylabel 'mach number: 𝓜    →'

plot for [i=1:words(datafiles)] word(datafiles,i) \
		u (scale_x($3)):(scale_y($4)) w l lw 2 t word(legends,i)
