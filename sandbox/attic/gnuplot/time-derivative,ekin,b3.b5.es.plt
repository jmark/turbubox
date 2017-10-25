set title 'Turbulent Box: Total Energy Turnover' offset 0,-0.5

set terminal svg enhanced
set encoding utf8

# note: you have to define crossing time t_c like 
#   $ gnuplot -e 't_c = ...' ...

crossingtime    = 5.8462062813e+13

  scaleX(x) = x / crossingtime 
 scaleX2(x) = x * crossingtime / (60*60*24*365*10**6)
scaleX2I(x) = x / crossingtime * (60*60*24*365*10**6)

  scaleY(x) = x / 10**33

set xtics 0.4
set xtics nomirror

# add second x axis for million years
set x2tics
set link x2 via scaleX2(x) inverse scaleX2I(x)

set xlabel  'crossing time {/:Italic t / t_c}  →'
set x2label 'time {/:Italic t} (Myr)   →'

set xrange [0:4]

set ylabel 'power rate {/:Italic P} ⋅ 10^{33} (erg/(ccm s))   →'

set grid

X = 1
Y = 3

plot \
    'data/cgs.128.b3,ekin.time-derivative' using (scaleX(column(X))):(scaleY(column(Y))) w l lw 2 t 'b3', \
    'data/cgs.128.b5,ekin.time-derivative' using (scaleX(column(X))):(scaleY(column(Y))) w l lw 2 t 'b5', \
    'data/cgs.128.es,ekin.time-derivative' using (scaleX(column(X))):(scaleY(column(Y))) w l lw 2 t 'es'
