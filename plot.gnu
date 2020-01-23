reset 
set terminal postscript enhanced 'Helvetica' 24 
set output 'plot.eps' 
plot [:][:] 'plot.dat' w p pt 13 ps 0.5 notitle 
