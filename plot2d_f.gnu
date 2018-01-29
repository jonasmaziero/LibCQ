reset 
set terminal postscript enhanced 'Helvetica' 24 
set output 'sigmoid.eps' 
 set title 'sigmoid function {/Symbol s}(x) = 1/(1+exp(-x))' 
set xlabel 'x' 
 set ylabel '{/Symbol s}' 
 plot [-10:10][:] 1/(1+exp(-x)) with lp notitle
