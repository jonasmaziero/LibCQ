reset 
set terminal postscript enhanced 'Helvetica' 24 
set output 'plot.eps' 
plot [:][0:1] 'plot.dat' u 1:2 w lp, '' u 1:3 w lp
