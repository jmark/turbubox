set terminal svg enhanced
set encoding utf8

# note: you have to define crossing time t_c like 
#   $ gnuplot -e 't_c = ...' ...

# cgs,mach=015
#crossingtime = 5.8462062813e+13

mach = 15
crossingtime = 0.06666667

# cgs,mach=100
#crossingtime = 8.7693094220e+12

  scaleX(x) = x / crossingtime 
 scaleX2(x) = x * crossingtime / (60*60*24*365*10**6)
scaleX2I(x) = x / crossingtime * (60*60*24*365*10**6)

  scaleY(x) = x / 10**46

set xtics 1
set xtics nomirror

# add second x axis for million years
set x2tics
set link x2 via scaleX2(x) inverse scaleX2I(x)

set xlabel  'crossing time {/:Italic t / t_c}  →'
set x2label 'time {/:Italic t} (Myr)   →'

set xrange [0:10]

set ylabel 'kinetic energy {/:Italic 𝓚} ⋅ 10^{46} (erg/ccm)   →'
set grid

plot \
    'data/unit_128_es_mach-015_evolution.dsv' using (scaleX($3)):(scaleY($4)) w l lw 2 t 'es'
