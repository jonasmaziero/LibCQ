//-----------------------------------------------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
//-----------------------------------------------------------------------------------------------------------------------------------
double f(double x, void *params){
  (void)(params); /* avoid unused parameter warning */
  return sin(x);
}
//-----------------------------------------------------------------------------------------------------------------------------------
void derivative_test(){
// compile with: gcc tests.c derivative.c gsl.c -lgsl -lgslcblas -lm
  gsl_function F;
  double result, abserr;
  F.function = &f;
  F.params = 0;
  
  FILE *fd = fopen("der.dat", "wr");
  double dx = 0.01, x = 0.0;
  while(x < 8.0*3.1415){
    x += dx;
    gsl_deriv_central (&F, x, 1e-8, &result, &abserr);
    fprintf(fd, "%f %f %f \n", x, result, cos(x));
  }
  
  // writing a gnuplot script, from C
  FILE *fg = fopen("der.gnu", "w");
  fprintf(fg, "reset \n");
  fprintf(fg, "set terminal postscript enhanced 'Helvetica' 24 \n");
  fprintf(fg, "set output 'der.eps' \n ");
  fprintf(fg, "plot 'der.dat' u 1:2 w p, '' u 1:3 w l \n");
  fclose(fg);
  // running a gnuplot script, from C
  system("gnuplot der.gnu");
  // seeing the result
  system("evince der.eps");
  //system("open -a skim der.eps");
}
//-----------------------------------------------------------------------------------------------------------------------------------
