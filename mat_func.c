#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

void proj(int *d, long double _Complex *psi, long double _Complex *pj) {
  int j, k;
  for (j = 0; j < (*d); j++) {
    for (k = 0; k < (*d); k++) {
      *(pj+j*(*d)+k) = (creall(*(psi+j))+I*cimagl(*(psi+j)))*(creall(*(psi+k))-I*cimagl(*(psi+k)));
    }
  }
}

long double _Complex inner_hs(int *d, long double _Complex *A, long double _Complex *B){
  int j,k;
  long double _Complex ip = 0;
  for(j = 0; j < (*d); j++){
    for(k = 0; k < (*d); k++){
      ip += (creall(*(A+j*(*d)+k))-I*cimagl(*(A+j*(*d)+k)))*(creall(*(B+j*(*d)+k))+I*cimagl(*(B+j*(*d)+k)));
    }
  }
  return ip;
}

long double norm_hs(int *d, long double _Complex *A){
  long double _Complex inner_hs(int *, long double _Complex *, long double _Complex *);
  long double norm = sqrtl(creall(inner_hs(d, A, A)));
  return norm;
}

long double trace(int *d, long double *A) {
  int j;
  long double tr = 0.0;
  for (j = 0; j < (*d); j++) {
    tr += *(A+j*(*d)+j);
  }
  return tr;
}

long double _Complex trace_c(int *d, long double _Complex *A) {
  int j;
  long double _Complex tr = 0.0;
  for (j = 0; j < (*d); j++) {
    tr += *(A+j*(*d)+j);
  }
  return tr;
}

void array_display(int *nr, int *nc, long double *A){
  int j,k;
  for(j = 0; j < (*nr); j++){
    for(k = 0; k < (*nc); k++){
      printf("%f \t", ((double) *(A+j*(*nc)+k)));
    }
    printf("\n");
  }
}

void array_display_c(int *nr, int *nc, long double _Complex *A){
  int j,k;
  printf("real part \n");
  for(j = 0; j < (*nr); j++){
    for(k = 0; k < (*nc); k++){
      printf("%f \t",((double) creall(*(A+j*(*nc)+k))));
    }
    printf("\n");
  }
  printf("imaginary part \n");
  for(j = 0; j < (*nr); j++){
    for(k = 0; k < (*nc); k++){
      printf("%f \t",((double) cimagl(*(A+j*(*nc)+k))));
    }
    printf("\n");
  }
}

long double veccsum(int *d, long double *vec) {
  long double vcs = 0.0;
  int j;
  for (j = 0; j < (*d); j++) {
    vcs += *(vec+j);
  }
  return vcs;
}

int veccsum_i(int *d, int *vec) {
  int vcs = 0;
  int j;
  for (j = 0; j < (*d); j++) {
    vcs += *(vec+j);
  }
  return vcs;
}

/*
void main() {
  long double _Complex *s1, *s2;
  s1 = (long double _Complex*)malloc(4*sizeof(long double _Complex));
  s2 = (long double _Complex*)malloc(4*sizeof(long double _Complex));
  *(s1+0*2+0) = 0.0; *(s1+0*2+1) = 1.0; *(s1+1*2+0) = 1.0; *(s1+1*2+1) = 0.0; 
  *(s2+0*2+0) = 0.0; *(s2+0*2+1) = -I; *(s2+1*2+0) = I; *(s2+1*2+1) = 0.0;
  void array_display_c(int *, int *, long double _Complex *);
  int nr = 2, nc = 2;
  array_display_c(&nr, &nc, s1); printf("\n");
  long double _Complex inner_hs(int *, long double _Complex *, long double _Complex *), ip;
  ip = inner_hs(&nc, s1, s1);
  printf("%Lf \n", creall(ip));
  long double norm_hs(int *, long double _Complex *), norm;
  norm = norm_hs(&nc, s1);
  printf("%Lf \n", norm);
}
*/

// gcc mat_func.c -lm