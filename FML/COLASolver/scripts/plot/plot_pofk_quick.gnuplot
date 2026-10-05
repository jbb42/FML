#!/usr/bin/env gnuplot

set term qt size 1800,1200 font ",12"
set logscale xy
set xlabel 'Wavenumber k [h/Mpc]'
set ylabel 'Power P_{ij}(k) [(Mpc/h)^3]'
set key top right box

dir = '/mn/stornext/u3/jonasbbe/pc/Dokumenter/FML/FML/COLASolver/output/snapshot_TestSim_z0.000/'

plot dir.'pofk_00.txt' using 1:2 with lines lw 2 title 'P_{00} (Matter Auto)', \
     dir.'pofk_01.txt' using 1:2 with lines lw 2 title 'P_{01} (Matter - \delta_L)', \
     dir.'pofk_11.txt' using 1:2 with lines lw 2 title 'P_{11} (\delta_L Auto)', \
     dir.'pofk_02.txt' using 1:2 with lines lw 2 title 'P_{02} (Matter - \delta_L^2)', \
     dir.'pofk_12.txt' using 1:2 with lines lw 2 title 'P_{12} (\delta_L - \delta_L^2)', \
     dir.'pofk_22.txt' using 1:2 with lines lw 2 title 'P_{22} (\delta_L^2 Auto)'
