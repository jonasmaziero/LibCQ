//-----------------------------------------------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
//-----------------------------------------------------------------------------------------------------------------------------------
void gnuplot_test(){
  // writing a gnuplot script, from C
  FILE *f = fopen("plot2d_f.gnu", "w");
  fprintf(f, "reset \n");
  fprintf(f, "set terminal postscript enhanced 'Helvetica' 24 \n");
  fprintf(f, "set output 'sigmoid.eps' \n ");
  fprintf(f, "set title 'sigmoid function {/Symbol s}(x) = 1/(1+exp(-x))' \n");
  fprintf(f, "set xlabel 'x' \n ");
  fprintf(f, "set ylabel '{/Symbol s}' \n ");
  fprintf(f, "plot [-10:10][:] 1/(1+exp(-x)) with lp notitle\n");
  fclose(f);
  // running a gnuplot script, from C
  system("gnuplot plot2d_f.gnu");
  // seeing the result
  system("evince sigmoid.eps");
}
//-----------------------------------------------------------------------------------------------------------------------------------
void plot2d(){ // plot data from a file
  double dx = 0.01, x = -10.0 - dx;
  double sigmoid(), sig;
  int j;
  FILE *fd = fopen("plot2d.dat", "wr");
  while (x < 10.0){
    x = x + dx;  sig = sigmoid(x);  fprintf(fd, "%f %f \n", x, sig);
  }
  // writing a gnuplot script, from C
  FILE *f = fopen("plot2d.gnu", "w");
  fprintf(f, "reset \n");
  fprintf(f, "set terminal postscript enhanced 'Helvetica' 24 \n");
  fprintf(f, "set output 'plot2d.eps' \n ");
  fprintf(f, "set xlabel 'x' \n ");
  fprintf(f, "set ylabel '{/Symbol s}' \n ");
  fprintf(f, "set xrange [-10:10] \n ");
  fprintf(f, "set yrange [-0.1:1.1] \n ");
  fprintf(f, "plot 'plot2d.dat' u 1:2 w p notitle, 1/(1+exp(-x)) w l notitle \n");
  fclose(f);
  // running a gnuplot script, from C
  system("gnuplot plot2d.gnu");
  // seeing the result
  system("evince plot2d.eps");
}
//-----------------------------------------------------------------------------------------------------------------------------------