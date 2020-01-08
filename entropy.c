
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

long double neumann(int *d, long double _Complex *rho){ 
  // Returns the von Neumann entropy of a density matrix
  double _Complex *A;
  A = (double _Complex *)malloc((*d)*(*d)*sizeof(double _Complex)); 
  int j, k;
  for (j = 0; j < (*d); j++) {
    for (k = 0; k < (*d); k++) {
      *(A+j*(*d)+k) = ((double _Complex) *(rho+j*(*d)+k));
    }
  }
  char jobz = 'N';
  double *W;
  W = (double *)malloc((*d)*sizeof(double)); 
  lapacke_zheevd(&jobz, d, A, W);
  long double vne, shannon(int *, long double *);
  long double *egva; egva = (long double *)malloc((*d)*sizeof(long double));
  for (j = 0; j < (*d); j++) {
    *(egva+j) = ((long double) *(W+j));
  }
  vne = shannon(d, egva);
  free(W); free(egva); free(A);
  return vne;
}

long double shannon(int *d, long double *pv) {
  int j;
  long double se = 0.0;
  for (j = 0; j < (*d); j++) { 
    if ((pv[j] > 1.e-15) && (pv[j] < (1.0-1.e-15))) { 
     se -= (*(pv+j))*log2l(*(pv+j));
    }  
  }
  return se;
}


/*
int main() {
  int d = 2;
  long double pv[] = {0.75,0.25};
  long double shannon(int *, long double *);  
  printf("%f \n", ((double) shannon(&d, pv)));

  long double _Complex *rho;
  rho = (long double _Complex *)malloc(d*d*sizeof(long double _Complex)); 
  *(rho+0*d+0) = 0.75; *(rho+0*d+1) = 0.0; *(rho+1*d+0) = 0.0; *(rho+1*d+1) = 0.25;
  long double neumann(int *, long double _Complex *);
  printf("%f \n", ((double) neumann(&d, rho)));
  return 0;
}
*/

// gcc entropy.c lapack.c -llapacke -lm