#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void rpv_zhsl(int *d, long double *rpv) {
  double genrand64_real1();
  long double norm;
  int j;
  *(rpv+0) = 1.0 - powl(((long double) genrand64_real1()),1.0/((*d)-1.0)); 
  norm = *(rpv+0);
  if ((*d) > 2) {
    for (j = 1; j < (*d)-1; j++) {
      *(rpv+j) = (1.0 - powl(((long double) genrand64_real1()),1.0/((*d)-j-1)))*(1.0-norm);
      norm = norm + (*(rpv+j));
    }
  }
  *(rpv+(*d)-1) = 1.0 - norm;
}

void rpv_devroye(int *d, long double *rpv) {
  void rng_exp(long double *);
  long double ern;
  int j;
  for (j = 0; j < (*d); j++) {
    rng_exp(&ern); *(rpv+j) = ern;
  }
  long double veccsum(int *, long double *), vcs; 
  vcs = veccsum(d, rpv);
  for (j = 0; j < (*d); j++) {
    *(rpv+j) /= vcs;
  } 
}

/*
int main() {
  void rng_init(); rng_init();
  int d = 10, ns = pow(10,4), ni = pow(10,2);
  long double delta = 1.0/((long double) ni);
  long double *rpva; rpva = (long double *)malloc(d*sizeof(long double));
  long double *rpv; rpv = (long double *)malloc(d*sizeof(long double));
  void rpv_zhsl(int *, long double *);
  void rpv_devroye(int *, long double *);
  int *ct; ct = (int *)malloc(ni*d*sizeof(int));
  int j,k,l;
  FILE *fd = fopen("plot.dat", "w");
  for (j = 0; j < ns; j++) {
    rpv_zhsl(&d, rpv); 
    //rpv_devroye(&d, rpv);
    for (k = 0; k < d; k++) {
      *(rpva+k) += *(rpv+k);
      for (l = 0; l < ni; l++) {   
        if (*(rpv+k) >= ((long double) l)*delta && *(rpv+k) < ((long double) (l+1))*delta) { 
          *(ct+l*d+k) += 1;
        } 
      }
    } 
  }
  if ( d <= 5 ) {
    for (j = 0; j < d; j++) {
      printf("%f \t", ((double) *(rpva+j))/((double) ns));
    }
    printf("\n");
  }
  for (l = 0; l < ni; l++) {
    fprintf(fd, "%f %f %f %f \n", ((double) l)*((double) delta), ((double) *(ct+l*d+1))/((double) ns), ((double) *(ct+l*d+2))/((double) ns), ((double) *(ct+l*d+3))/((double) ns));
  }
  fclose(fd);
  void plot(); plot();
  free(rpv); free(rpva); free(ct);
  return 0;
}
*/


// gcc rpvg.c rng.c MT19937_64.c mat_func.c plots.c -lm

