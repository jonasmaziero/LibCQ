#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void rpv_zhsl(int *d, double *rpv) {
  double genrand64_real1(), norm;
  int j;
  *(rpv+0) = 1.0 - pow(genrand64_real1(),1.0/((*d)-1.0)); norm = *(rpv+0);
  if ((*d) > 2) {
    for (j = 1; j < (*d)-1; j++) {
      *(rpv+j) = (1.0 - pow(genrand64_real1(),1.0/((*d)-j-1)))*(1.0-norm);
      norm = norm + (*(rpv+j));
    }
  }
  *(rpv+(*d)-1) = 1.0 - norm;
}

void rpv_devroye(int *d, double *rpv) {
  void rng_exp(double *);
  double ern;
  int j;
  for (j = 0; j < (*d); j++) {
    rng_exp(&ern); *(rpv+j) = ern;
  }
  double veccsum(int *, double *), vcs; 
  vcs = veccsum(d, rpv);
  for (j = 0; j < (*d); j++) {
    *(rpv+j) /= vcs;
  } 
}

/*
int main() {
  unsigned long long int seed;  srand(time(NULL));
  seed = rand();  init_genrand64(seed);
  int d = 100, ns = pow(10,4), ni = pow(10,2);
  double delta = 1.0/((double) ni);
  double *rpva; rpva = (double *)malloc(d*sizeof(double));
  double *rpv; rpv = (double *)malloc(d*sizeof(double));
  void rpv_zhsl(int *, double *);
  void rpv_devroye(int *, double *);
  int *ct; ct = (int *)malloc(ni*d*sizeof(int));
  int j,k,l;
  FILE *fd = fopen("plot.dat", "w");
  for (j = 0; j < ns; j++) {
    //rpv_zhsl(&d, rpv); 
    rpv_devroye(&d, rpv);
    for (k = 0; k < d; k++) {
      *(rpva+k) += *(rpv+k);
      //if (*(rpv+k) == 1.0) *(rpv+k) -= 1.e-10;
      for (l = 0; l < ni; l++) {   
        if (*(rpv+k) >= ((double) l)*delta && *(rpv+k) < ((double) (l+1))*delta) { 
          *(ct+l*d+k) += 1;
        } 
      }
    } 
  }
  if ( d <= 5 ) {
    for (j = 0; j < d; j++) {
      printf("%f \t", *(rpva+j)/((double) ns));
    }
    printf("\n");
  }
  for (l = 0; l < ni; l++) {
    fprintf(fd, "%f %f %f %f \n", ((double) l)*delta, ((double) *(ct+l*d+1))/((double) ns), ((double) *(ct+l*d+2))/((double) ns), ((double) *(ct+l*d+3))/((double) ns));
  }
  fclose(fd);
  FILE *fg = fopen("plot.gnu", "w");
  fprintf(fg, "reset \n");
  fprintf(fg, "set terminal postscript enhanced 'Helvetica' 24 \n");
  fprintf(fg, "set output 'plot.eps' \n ");
  fprintf(fg, "plot [0:1][0:] 'plot.dat' u 1:2 w lp,'' u 1:3 w lp,'' u 1:4 w lp \n");
  fclose(fg);
  system("gnuplot plot.gnu");
  system("evince plot.eps & \n");
  free(rpv); free(rpva); free(ct);
  return 0;
}
*/

// gcc rpvg.c rng.c MT19937_64.c mat_func.c -lm

