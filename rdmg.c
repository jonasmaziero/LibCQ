// gcc main.c rdmg.c MT19937_64.c -lm
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "MT19937_64.h"
#include <sys/utsname.h>

void ginibre(int *d, double _Complex *G){

 G = np.zeros((d, d), dtype=complex)
    mu, sigma = 0.0, 1.0
    for j in range(0, d):
        grn = np.random.normal(mu, sigma, 2*d)
        for k in range(0, d):
            G[j][k] = grn[k] + (1j)*grn[k+d]
}

void test_rand_circle(){
  int n, j;
  double r, x, y;
  void rand_circle(double *,double *,double *);
  FILE *fd = fopen("plot.dat", "w");
  n = 100000; r = 1.0;
  for (j = 0; j < n; j++){
    rand_circle(&r,&x,&y);
    fprintf(fd, "%f %f \n", x, y);
  }
  fclose(fd);
  FILE *fg = fopen("plot.gnu", "w");
  fprintf(fg, "reset \n");
  fprintf(fg, "set terminal postscript enhanced 'Helvetica' 24 \n");
  fprintf(fg, "set output 'plot.eps' \n ");
  fprintf(fg, "plot 'plot.dat' w p pt 7 ps 0.2 lc 3 \n");
  fclose(fg);
  system("gnuplot plot.gnu");
  struct utsname osf;  uname(&osf);
  long unsigned int strlen();
  if(strlen(osf.sysname) == 5){
    system("evince plot.eps & \n");
  } else if(strlen(osf.sysname) == 6){
    system("open -a skim plot.eps & \n");
  }
}

void rand_circle(double *r, double *x, double *y){
  double th, ph, rh, pi = M_PI;
  double genrand64_real1();
  th = 2*pi*genrand64_real1();
  rh = (*r)*sqrt(genrand64_real1());
  (*x) = rh*cos(th);
  (*y) = rh*sin(th);
}
