// gcc entropy.c lapack.c -llapacke -lm
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double neumann(int *d, double _Complex *rho){ 
  // Returns the von Neumann entropy of a density matrix
  double _Complex *A;
  A = (double _Complex *)malloc((*d)*(*d)*sizeof(double _Complex)); 
  int j, k;
  for (j = 0; j < (*d); j++) {
    for (k = 0; k < (*d); k++) {
      *(A+j*(*d)+k) = *(rho+j*(*d)+k);
    }
  }
  char jobz = 'N';
  double *W;
  W = (double *)malloc((*d)*sizeof(double)); 
  lapacke_zheevd(&jobz, d, A, W);
  free(A);
  double vne, shannon(int *, double *);
  vne = shannon(d, W);
  free(W);
  return vne;
}

double shannon(int *d, double *pv) {
  int j;
  double se = 0.0;
  for (j = 0; j < (*d); j++) { 
    if ((pv[j] > 1.e-15) && (pv[j] < (1.0-1.e-15))) { 
     se -= (*(pv+j))*log2(*(pv+j));
    }  
  }
  return se;
}

/*
int main() {
  int d = 2;
  double pv[] = {0.5,0.5};
  double shannon(int *, double *);  
  printf("%f \n", shannon(&d, pv));

  double _Complex *rho;
  rho = (double _Complex *)malloc(d*d*sizeof(double _Complex)); 
  *(rho+0*d+0) = 0.75; *(rho+0*d+1) = 0.0; *(rho+1*d+0) = 0.0; *(rho+1*d+1) = 0.25;
  double neumann(int *, double _Complex *);
  printf("%f \n", neumann(&d, rho));
  return 0;
}
*/