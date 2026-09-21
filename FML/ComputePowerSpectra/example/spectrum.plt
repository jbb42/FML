set term qt size 1800,1200 font ",12"
set logscale xy
set xlabel 'Wavenumber k'
set ylabel 'Power P(k)'
set key top right box

plot 'my_spectrum.txt' every :::0::0 using 1:2 with linespoints pt 7 lw 2 title 'Standard P(k)', \
     'my_spectrum.txt' every :::1::1 using 1:2 with linespoints pt 5 lw 2 title 'Interlaced P(k)'
