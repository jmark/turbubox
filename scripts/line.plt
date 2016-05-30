set grid

set style line 1 lw 2 pt 3 ps 0.5
set style line 2 lw 2 pt 3 ps 0.5
set style line 3 lw 2 pt 3 ps 0.5



set title "Line slice: bouchut5"

plot \
    'line-z.csv'    using 6:7 w l t "dens" ls 1, \
    'line-z.csv'    using 6:8 w l t "velx" ls 2
 
