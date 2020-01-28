#include <stdio.h>
#include <stdlib.h>

void plot() {
  FILE *fgnu = fopen("plot.gnu", "w");
  fprintf(fgnu, "reset \n");
  fprintf(fgnu, "set terminal postscript enhanced 'Helvetica' 24 \n");
  fprintf(fgnu, "set output 'plot.eps' \n");
  //fprintf(fgnu, "plot [:][:] 'plot.dat' w p pt 13 ps 0.5 notitle \n");
  //fprintf(fgnu, "plot [:][0:1] 'plot.dat' u 1:2 w lp, '' u 1:3 w lp, '' u 1:4 w lp\n");
  fprintf(fgnu, "plot [:][0:1] 'plot.dat' u 1:2 w lp, '' u 1:3 w lp\n");
  //fprintf(fgnu, "plot 'plot.dat' w p pt 13 ps 0.5 notitle, exp(-x) \n");
  //fprintf(fgnu, "plot 'plot.dat' w p notitle, (1/sqrt(2*pi))*exp(-0.5*x**2)\n");
  fclose(fgnu);
  system("gnuplot plot.gnu");
  system("evince plot.eps&");
}