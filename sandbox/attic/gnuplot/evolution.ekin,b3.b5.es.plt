set title 'Turbulent Box: Total Kinetic Energy'

set terminal svg enhanced
set encoding utf8

# note: you have to define crossing time t_c like 
#   $ gnuplot -e 't_c = ...' ...

crossingtime    = 5.8462062813e+13

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

set xrange [0:6]

set ylabel 'kinetic energy {/:Italic 𝓚} ⋅ 10^{46} (erg/ccm)   →'

set grid

plot \
    'data/cgs.128.b3,evolution' using (scaleX($3)):(scaleY($4)) w l lw 1 t 'b3', \
    'data/cgs.128.b5,evolution' using (scaleX($3)):(scaleY($4)) w l lw 1 t 'b5', \
    'data/cgs.128.es,alpha=1.0,flux=2,evolution' using (scaleX($3)):(scaleY($4)) w l lw 1 t 'es: α=1.000, flux=2', \
    'data/cgs.128.es,alpha=1.000,evolution' using (scaleX($3)):(scaleY($4)) w l lw 1 t 'es: α=1.000, flux=3'
