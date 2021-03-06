#set title 'Turbulent Box: Evolution of the Mach Number'

set terminal svg enhanced
set encoding utf8

# note: you have to define crossing time t_c like 
#   $ gnuplot -e 't_c = ...' ...

# Following values refer to Mach-15 stirring:

# cgs,mach=015
#crossingtime = 5.8462062813e+13

# cgs,mach=100
#crossingtime = 8.7693094220e+12

#domainsize      = 3.0857000000e+19
#ambientdens     = 1.6726218980e-23
#soundspeed      = 3.5187491415e+04

# unit,mach=100
unitsystem='unit'
gridsize='128'
solver='b3'
mach='015'

crossingtime = 0.06666667
domainsize      = 1.0
ambientdens     = 1.0
soundspeed      = 1.0

totalmass       = domainsize**3 * ambientdens

#  scaleX(x) = x / crossingtime 
# scaleX2(x) = x * crossingtime / (60*60*24*365*10**6)
#scaleX2I(x) = x / crossingtime * (60*60*24*365*10**6)

  scaleX(x) = x / crossingtime
 scaleX2(x) = x
scaleX2I(x) = x

  scaleY(x) = sqrt(2.0 * x / totalmass) / soundspeed


set xtics 1

set xtics nomirror

# add second x axis for million years
#set x2tics
#set link x2 via scaleX2(x) inverse scaleX2I(x)

set xlabel  'crossing time: {/:Italic t / t_c}   →'
#set x2label 'time: {/:Italic t} (Myr)   →'

#set xrange [0:10]

set ytics 1
set ylabel 'mach number: 𝓜    →'

set grid


plot \
    'data/unit_128_es_mach-'.mach.'_evolution.dsv' using (scaleX($3)):(scaleY($4)) w l lw 2 t 'es'

    #'data/unit,128,b3,mach='.mach.',evolution.dsv' using (scaleX($3)):(scaleY($4)) w l lw 2 t 'b3',\
    #'data/unit,128,b5,mach='.mach.',evolution.dsv' using (scaleX($3)):(scaleY($4)) w l lw 2 t 'b5',\
    #'data/unit,128,es,mach='.mach.',evolution.dsv' using (scaleX($3)):(scaleY($4)) w l lw 2 t 'es'
    #'data/cgs,128,b5,mach='.mach.',evolution.dsv' using (scaleX($3)):(scaleY($4)) w l lw 2 t 'b5',\
    #'data/cgs,128,es,mach='.mach.',evolution.dsv' using (scaleX($3)):(scaleY($4)) w l lw 2 t 'es'
