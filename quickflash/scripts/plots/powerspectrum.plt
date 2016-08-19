set grid

set style line 1 lw 2 pt 3 ps 0.5
set style line 2 lw 2 pt 3 ps 0.5
set style line 3 lw 2 pt 3 ps 0.5

set title "Powerspectrum"

set xlabel "k ->"
set x2label "log10(k) ->"
set ylabel "log10(P(k)dk) ->"

set logscale xy 10

set x2tics nomirror

set xtics 0,0.5,70

set yrange [0.25:3.1]

plot \
    '0099_8waveclean'    w l ls 1 t '8wave', \
    '0099_bouchut5clean' w l ls 3 t 'bouchut5',\
    '0099_bouchut3clean' w l ls 2 t 'bouchut3'
