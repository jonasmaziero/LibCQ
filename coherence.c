#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

long double coh_l1(int *d, long double _Complex *rho) {  
  // l1-norm coherence, Ref: PRL 113, 140401 (2014)
  int j, k;
  long double coh = 0.0;
  for (j = 0; j < ((*d)-1); j++) {
    for (k = j+1; k < (*d); k++) { 
      coh += sqrtl(powl(creall(*(rho+j*(*d)+k)),2.0) + powl(cimagl(*(rho+j*(*d)+k)),2.0));
    }
  }
  return 2.0*coh/((long double) ((*d)-1));
}

long double coh_hs(int *d, long double _Complex *rho) {
  int j, k;
  long double coh = 0.0;
  for (j = 0; j < ((*d)-1); j++) {
    for (k = j+1; k < (*d); k++) {  
      coh += powl(creall(*(rho+j*(*d)+k)),2.0) + powl(cimagl(*(rho+j*(*d)+k)),2.0);
    }
  }
  return 2.0*coh/(((long double) ((*d)-1))/((long double) (*d)));
}

long double coh_re(int *d, long double _Complex *rho) { 
  // relative entropy of quantum coherence, Ref: PRL 113, 140401 (2014)
  int j;
  long double *pv;
  pv = (long double *)malloc((*d)*sizeof(long double));
  for(j = 0; j < (*d); j++){
    *(pv+j) = creall(*(rho+j*(*d)+j));
  }
  long double coh, shannon(int *, long double *);
  long double neumann(int *, long double _Complex *);
  coh = shannon(d, pv) - neumann(d, rho);
  free(pv);
  return coh/log2l((*d));
}

/*
int main() {
  long double theta = 0.0, phi = 0.0, delta = 0.05;
  int d = 2;
  long double _Complex *psi;
  psi = (long double _Complex *)malloc(d*sizeof(long double _Complex));
  long double _Complex *rho;
  rho = (long double _Complex *)malloc(d*d*sizeof(long double _Complex));
  long double coh1, coh2;
  long double coh_l1(int *, long double _Complex *);
  long double coh_re(int *, long double _Complex *);
  void psi1qb(long double *, long double *, long double _Complex *);
  void proj(int *, long double _Complex *, long double _Complex *);
  theta = 0.0;
  while (theta < M_PI) {
    psi1qb(&theta, &phi, psi);  
    proj(&d, psi, rho); 
    coh1 = coh_l1(&d, rho);  coh2 = coh_re(&d, rho);
    printf("%f \t %f \t %f \n", ((double) theta), ((double) coh1), ((double) coh2));
    theta += delta;
  }
  free(psi);
  free(rho);
  return 0;
}
*/

// gcc coherence.c entropy.c lapack.c states.c mat_func.c -llapacke -lm