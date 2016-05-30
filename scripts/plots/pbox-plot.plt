set terminal pngcairo size 500,500

set title 'MHD Gir. Turb. 64x64 at z = 0, n = $NR'

#set xlabel 'x ->'
#set ylabel 'y ->'

set lmargin at screen 0
set rmargin at screen 1
set tmargin at screen 1
set bmargin at screen 0
set xtics offset 0, screen 0.1
set ytics offset screen 0.1, 0

unset colorbox
unset xlabel
unset ylabel

set size ratio 1
set view map

splot '$FILE' u 3:4:5 with pm3d t 'n = $NR'
