set xlabel "Position x"
set ylabel "Deviation"
plot "plot.dat" u 1:2 title 'Solved' w l
replot -9.0/8+cosh(x-0.5) title 'Expected'
pause -1
