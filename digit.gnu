reset 
set terminal postscript enhanced 'Helvetica' 24 
set output 'digit.eps' 
set pm3d map clip4in corners2col c1 
splot 'digit.dat' matrix using 1:(1-$2):3 with pm3d notitle 
