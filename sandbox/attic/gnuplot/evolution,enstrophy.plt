set terminal svg enhanced
set encoding utf8

# note: you have to define crossing time t_c like 
#   $ gnuplot -e 't_c = ...' ...

# cgs,mach=015
#crossingtime = 5.8462062813e+13

# cgs,mach=100
#crossingtime = 8.7693094220e+12

# unit,mach=015
crossingtime = 0.06666667

  scaleX(x) = x / crossingtime 
 scaleX2(x) = x * crossingtime / (60*60*24*365*10**6)
scaleX2I(x) = x / crossingtime * (60*60*24*365*10**6)

  scaleY(x) = x / 10**44

set xtics 1
set xtics nomirror

# add second x axis for million years
set x2tics
set link x2 via scaleX2(x) inverse scaleX2I(x)

set xlabel  'crossing time {/:Italic t / t_c}  →'
set x2label 'time {/:Italic t} (Myr)   →'

set xrange [0:10]

set ylabel 'enstrophy {/:Italic 𝓔} ⋅ 10^{44} (erg/ccm)   →'

set grid

mach='100'

plot \
    'data/cgs,128,b3,mach='.mach.',evolution.dsv' using (scaleX($3)):(scaleY($5)) w l lw 2 t 'b3',\
    'data/cgs,128,b5,mach='.mach.',evolution.dsv' using (scaleX($3)):(scaleY($5)) w l lw 2 t 'b5',\
    'data/cgs,128,es,mach='.mach.',evolution.dsv' using (scaleX($3)):(scaleY($5)) w l lw 2 t 'es'
