#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

double coh_l1(int *d, double _Complex *rho) {
  // l1-norm coherence, Ref: PRL 113, 140401 (2014)
  int j, k;
  double coh = 0.0;
  for (j = 0; j < ((*d)-1); j++) {
    for (k = j+1; k < (*d); k++) { 
      coh += sqrt(pow(creal(*(rho+j*(*d)+k)),2.0) + pow(cimag(*(rho+j*(*d)+k)),2.0));
    }
  }
  return 2.0*coh/((double) ((*d)-1));
}

double coh_hs(int *d, double _Complex *rho) { // normalized
  int j, k;
  double coh = 0.0;
  for (j = 0; j < ((*d)-1); j++) {
    for (k = j+1; k < (*d); k++) {  
      coh += pow(creal(*(rho+j*(*d)+k)),2) + pow(cimag(*(rho+j*(*d)+k)),2);
    }
  }
  coh = 2*coh/(((double) ((*d)-1))/((double) (*d)));
  return coh;
}

double coh_hs_(int *d, double _Complex *rho) { // non-normalized
  int j, k;
  double coh = 0.0;
  for (j = 0; j < ((*d)-1); j++) {
    for (k = j+1; k < (*d); k++) {  
      coh += pow(creal(*(rho+j*(*d)+k)),2) + pow(cimag(*(rho+j*(*d)+k)),2);
    }
  }
  return 2*coh;
}

double coh_re(int *d, double _Complex *rho) { 
  // relative entropy of quantum coherence, Ref: PRL 113, 140401 (2014)
  int j;
  double *pv;
  pv = (double *)malloc((*d)*sizeof(double));
  for(j = 0; j < (*d); j++){
    *(pv+j) = creal(*(rho+j*(*d)+j));
  }
  double coh, shannon(int *, double *);
  double neumann(int *, double _Complex *);
  coh = shannon(d, pv) - neumann(d, rho);
  free(pv);
  return coh/log2((*d));
}

/*
int main() {
  double theta = 0.0, phi = 0.0, delta = 0.05;
  int d = 2;
  double _Complex *psi;
  psi = (double _Complex *)malloc(d*sizeof(double _Complex));
  double _Complex *rho;
  rho = (double _Complex *)malloc(d*d*sizeof(double _Complex));
  double coh1, coh2, coh3;
  double coh_l1(int *, double _Complex *);
  double coh_hs(int *, double _Complex *);
  double coh_re(int *, double _Complex *);
  void psi1qb(double *, double *, double _Complex *);
  void proj(int *, double _Complex *, double _Complex *);
  FILE *fd = fopen("plot.dat", "w");
  theta = 0.0;
  while (theta < M_PI) {
    psi1qb(&theta, &phi, psi);  
    proj(&d, psi, rho); 
    coh1 = coh_l1(&d, rho); coh2 = coh_hs(&d, rho); coh3 = coh_re(&d, rho);
    fprintf(fd,"%f %f %f %f \n", ((double) theta), ((double) coh1), ((double) coh2), ((double) coh3));
    theta += delta;
  }
  free(psi); free(rho);
  fclose(fd);
  void plot(); plot();
  return 0;
}
*/

// gcc coherence.c entropy.c lapack.c states.c mat_func.c plots.c -llapacke -lm