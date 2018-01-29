reset 
set terminal postscript enhanced 'Helvetica' 24 
set output 'rn.eps' 
plot 'rn.dat', (1.0/sqrt(2.0*pi*(1.0**2)))*exp(-((x-0)**2)/(2.0*(1.0**2))) 
