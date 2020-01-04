#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

void proj(int *d, double _Complex *psi, double _Complex *pj) {
  int j, k;
  for (j = 0; j < (*d); j++) {
    for (k = 0; k < (*d); k++) {
      *(pj+j*(*d)+k) = (creal(*(psi+j))+I*cimag(*(psi+j)))*(creal(*(psi+k))-I*cimag(*(psi+k)));
    }
  }
}

double _Complex inner_hs(int *d, double _Complex *A, double _Complex *B){
  int j,k;
  double _Complex ip = 0;
  for(j = 0; j < (*d); j++){
    for(k = 0; k < (*d); k++){
      ip += (creal(*(A+j*(*d)+k))-I*cimag(*(A+j*(*d)+k)))*(creal(*(B+j*(*d)+k))+I*cimag(*(B+j*(*d)+k)));
    }
  }
  return ip;
}

double norm_hs(int *d, double _Complex *A){
  double _Complex inner_hs(int *, double _Complex *, double _Complex *);
  double norm = sqrt(creal(inner_hs(d, A, A)));
  return norm;
}

double trace(int *d, double *A) {
  int j;
  double tr = 0.0;
  for (j = 0; j < (*d); j++) {
    tr += *(A+j*(*d)+j);
  }
  return tr;
}

double _Complex trace_c(int *d, double _Complex *A) {
  int j;
  double _Complex tr = 0.0;
  for (j = 0; j < (*d); j++) {
    tr += *(A+j*(*d)+j);
  }
  return tr;
}

void array_display(int *nr, int *nc, double *A){
  int j,k;
  for(j = 0; j < (*nr); j++){
    for(k = 0; k < (*nc); k++){
      printf("%f \t",*(A+j*(*nc)+k));
    }
    printf("\n");
  }
}

void array_display_c(int *nr, int *nc, double _Complex *A){
  int j,k;
  printf("real part \n");
  for(j = 0; j < (*nr); j++){
    for(k = 0; k < (*nc); k++){
      printf("%f \t",creal(*(A+j*(*nc)+k)));
    }
    printf("\n");
  }
  printf("imaginary part \n");
  for(j = 0; j < (*nr); j++){
    for(k = 0; k < (*nc); k++){
      printf("%f \t",cimag(*(A+j*(*nc)+k)));
    }
    printf("\n");
  }
}

double veccsum(int *d, double *vec) {
  double vcs = 0.0;
  int j;
  for (j = 0; j < (*d); j++) {
    vcs += *(vec+j);
  }
  return vcs;
}

/*
void test_norm() {
  double _Complex *s1, *s2;
  s1 = (double _Complex*)malloc(4*sizeof(double _Complex));
  s2 = (double _Complex*)malloc(4*sizeof(double _Complex));
  *(s1+0*2+0) = 0.0; *(s1+0*2+1) = 1.0; *(s1+1*2+0) = 1.0; *(s1+1*2+1) = 0.0; 
  *(s2+0*2+0) = 0.0; *(s2+0*2+1) = -I; *(s2+1*2+0) = I; *(s2+1*2+1) = 0.0;
  void array_display_c(int *, int *, double _Complex *);
  int nr = 2, nc = 2;
  //array_display_c(&nr, &nc, s1);
  double _Complex inner_hs(int *, double _Complex *, double _Complex *), ip;
  ip = inner_hs(&nc, s1, s1);
  printf("%f \n", creal(ip));
  double norm_hs(int *, double _Complex *), norm;
  norm = norm_hs(&nc, s1);
  printf("%f \n", norm);
}
*/