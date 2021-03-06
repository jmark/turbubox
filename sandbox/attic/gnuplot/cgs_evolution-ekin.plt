#!/bin/env gnuplot

set terminal svg enhanced
set encoding utf8

datafiles	= "`echo $GPI_datafiles`"
legends		= "`echo $GPI_legends`"

mach		= `echo $GPI_mach`
crosstime	= `echo $GPI_crosstime`

scale_x(x)		= x / crosstime 
scale_x2(x)		= x * crosstime / (60*60*24*365*10**6)
scale_x2I(x)	= x / crosstime * (60*60*24*365*10**6)

scale_y(x)		= x / 10**46

set grid

set xtics nomirror
set x2tics
set link x2 via scale_x2(x) inverse scale_x2I(x)

set xrange [0:5]

set  xlabel 'crossing time {/:Italic t / t_c}  →'
set x2label 'time {/:Italic t} (Myr)   →'
set ylabel  'kinetic energy {/:Italic 𝓚} ⋅ 10^{46} (erg/ccm)   →'

plot for [i=1:words(datafiles)] word(datafiles,i) \
		u (scale_x($3)):(scale_y($4)) w l lw 2 t word(legends,i)
