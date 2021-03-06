set title 'Turbulent Box: Evolution of the Mach Number'

set terminal svg enhanced
set encoding utf8

# note: you have to define crossing time t_c like 
#   $ gnuplot -e 't_c = ...' ...

# Following values refer to Mach-15 stirring:

crossingtime    = 5.8462062813e+13
domainsize      = 3.0857000000e+19
ambientdens     = 1.6726218980e-23
soundspeed      = 3.5187491415e+04
totalmass       = domainsize**3 * ambientdens

  scaleX(x) = x / crossingtime 
 scaleX2(x) = x * crossingtime / (60*60*24*365*10**6)
scaleX2I(x) = x / crossingtime * (60*60*24*365*10**6)

  scaleY(x) = sqrt(2.0 * x / totalmass) / soundspeed

set xtics 1
set xtics nomirror

# add second x axis for million years
set x2tics
set link x2 via scaleX2(x) inverse scaleX2I(x)

set xlabel  'crossing time: {/:Italic t / t_c}   →'
set x2label 'time: {/:Italic t} (Myr)   →'

set xrange [0:6]

set ytics 1
set ylabel 'mach number: 𝓜    →'

set grid

plot \
    'data/cgs.128.es,alpha=0.500,evolution' using (scaleX($3)):(scaleY($4)) w l lw 1 t 'α=0.500, flux=3', \
    'data/cgs.128.es,alpha=0.750,evolution' using (scaleX($3)):(scaleY($4)) w l lw 1 t 'α=0.750, flux=3', \
    'data/cgs.128.es,alpha=0.875,evolution' using (scaleX($3)):(scaleY($4)) w l lw 1 t 'α=0.875, flux=3', \
    'data/cgs.128.es,alpha=1.000,evolution' using (scaleX($3)):(scaleY($4)) w l lw 1 t 'α=1.000, flux=3', \
    'data/cgs.128.es,alpha=1.0,flux=2,evolution' using (scaleX($3)):(scaleY($4)) w l lw 1 t 'α=1.000, flux=2'
