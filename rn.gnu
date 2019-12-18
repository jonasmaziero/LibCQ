reset 
set terminal postscript enhanced 'Helvetica' 24 
set output 'rn.eps' 
plot 'rn.dat' t 'numerico',(1/sqrt(2*pi))*exp(-x**2/2)  t 'exato' 
