
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

//#define I _Complex_I;


void proj(int *d, double _Complex *psi, double _Complex *pj) {
  int j, k;
  for (j = 0; j < (*d); j++) {
    for (k = 0; k < (*d); k++) {
      *(pj+j*(*d)+k) = (creal(*(psi+j))+I*cimag(*(psi+j)))*(creal(*(psi+k))-I*cimag(*(psi+k)));
    }
  }
}


double ip(int *d, double *v, double *w) {
  double ipd = 0.0;
  int j;
  for (j = 0; j < (*d); j++) {
    ipd += (*(v+j))*(*(w+j));
  }
  return ipd;
}


double _Complex ip_c(int *d, double _Complex *v, double _Complex *w) {
  double _Complex ip = 0;
  int j;
  for (j = 0; j < (*d); j++) {
    ip += conj(*(v+j))*(*(w+j));
  }
  return ip;
}

long double _Complex ip_cl(int *d, long double _Complex *v, long double _Complex *w) {
  long double _Complex ip = 0;
  int j;
  for (j = 0; j < (*d); j++) {
    ip += conjl(*(v+j))*(*(w+j));
  }
  return ip;
}


double norm(int *d, double *v) {
  double ip(int *, double *, double *);
  double n = sqrt(ip(d, v, v));
  return n;
}


double norm_c(int *d, double _Complex *v) {
  double _Complex ip_c(int *, double _Complex *, double _Complex *);
  double _Complex ipc = ip_c(d, v, v);
  double n = sqrt(pow(creal(ipc),2)+pow(cimag(ipc),2)); 
  return n;
}

double vnorm2_c(int *d, double _Complex *v) {
  int j;
  double vn = 0;
  for (j = 0; j < (*d); j++) {
    vn += (pow(creal(*(v+j)),2)+pow(cimag(*(v+j)),2));
  }
  return vn;
}


long double vnorm2_cl(int *d, long double _Complex *v) {
  int j;
  long double vn = 0;
  for (j = 0; j < (*d); j++) {
    vn += (powl(creall(*(v+j)),2)+powl(cimagl(*(v+j)),2));
  }
  return vn;
}


double _Complex ip_hs(int *d, double _Complex *A, double _Complex *B) {
  int j,k;
  double _Complex ip = 0;
  for (j = 0; j < (*d); j++) {
    for (k = 0; k < (*d); k++) {
      ip += (creal(*(A+j*(*d)+k))*creal(*(B+j*(*d)+k)));
      ip += (cimag(*(A+j*(*d)+k))*cimag(*(B+j*(*d)+k)));
      ip += I*(creal(*(A+j*(*d)+k))*cimag(*(B+j*(*d)+k)));
      ip -= I*(cimag(*(A+j*(*d)+k))*creal(*(B+j*(*d)+k)));
    }
  }
  return ip;
}


double norm_hs(int *d, double _Complex *A) {
  int j, k;
  double hsn = 0;
  for (j = 0; j < (*d); j++) {
    for (k = 0; k < (*d); k++) {
      hsn += (pow(creal(*(A+j*(*d)+k)),2) + pow(cimag(*(A+j*(*d)+k)),2));
    }
  }
  return sqrt(hsn);
}

double norm_hsl(int *d, long double _Complex *A) {
  int j, k;
  long double hsn = 0;
  for (j = 0; j < (*d); j++) {
    for (k = 0; k < (*d); k++) {
      hsn += (powl(creall(*(A+j*(*d)+k)),2) + powl(cimagl(*(A+j*(*d)+k)),2));
    }
  }
  return sqrtl(hsn);
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


void array_display(int *nr, int *nc, double *A) {
  int j,k;
  for(j = 0; j < (*nr); j++){
    for(k = 0; k < (*nc); k++){
      printf("%f \t", ((double) *(A+j*(*nc)+k)));
    }
    printf("\n");
  }
}


void array_display_c(int *nr, int *nc, double _Complex *A) {
  int j,k;
  printf("real part \n");
  for(j = 0; j < (*nr); j++){
    for(k = 0; k < (*nc); k++){
      printf("%.3f \t",((double) creal(*(A+j*(*nc)+k))));
    }
    printf("\n");
  }
  printf("imaginary part \n");
  for(j = 0; j < (*nr); j++){
    for(k = 0; k < (*nc); k++){
      printf("%.3f \t",((double) cimag(*(A+j*(*nc)+k))));
    }
    printf("\n");
  }
}

void array_display_i(int *nr, int *nc, int *A) {
  int j,k;
  for(j = 0; j < (*nr); j++){
    for(k = 0; k < (*nc); k++){
      printf("%i \t", *(A+j*(*nc)+k));
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


double veccsuml(int *d, long double *vec) { // long double version
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


void zero_mat(int *nr, int *nc, double *A) {
  int j, k;
  for (j = 0; j < (*nr); j++) {
    for (k = 0; k < (*nc); k++) {
      *(A+j*(*nc)+k) = 0;
    }
  }
}


void zero_mat_c(int *nr, int *nc, double _Complex *A) {
  int j, k;
  for (j = 0; j < (*nr); j++) {
    for (k = 0; k < (*nc); k++) {
      *(A+j*(*nc)+k) = 0;
    }
  }
}


void zero_mat_i(int *nr, int *nc, int *A) {
  int j, k;
  for (j = 0; j < (*nr); j++) {
    for (k = 0; k < (*nc); k++) {
      *(A+j*(*nc)+k) = 0;
    }
  }
}

void one_mat_i(int *nr, int *nc, int *A) {
  int j, k;
  for (j = 0; j < (*nr); j++) {
    for (k = 0; k < (*nc); k++) {
      *(A+j*(*nc)+k) = 1;
    }
  }
}


/*
int main() {
  double _Complex *s1, *s2;
  s1 = (double _Complex*)malloc(4*sizeof(double _Complex));
  s2 = (double _Complex*)malloc(4*sizeof(double _Complex));
  *(s1+0*2+0) = 0.0; *(s1+0*2+1) = 1.0; *(s1+1*2+0) = 1.0; *(s1+1*2+1) = 0.0; 
  *(s2+0*2+0) = 0.0; *(s2+0*2+1) = -I; *(s2+1*2+0) = I; *(s2+1*2+1) = 0.0;
  void array_display_c(int *, int *, double _Complex *);
  int nr = 2, nc = 2;
  array_display_c(&nr, &nc, s1); printf("\n");
  double _Complex inner_hs(int *, double _Complex *, double _Complex *), ip;
  ip = inner_hs(&nc, s1, s1);
  printf("%f \n", creal(ip));
  double norm_hs(int *, double _Complex *), norm;
  norm = norm_hs(&nc, s1);
  printf("%f \n", norm);
  int d = 2;
  double a[d], b[d]; a[0] = 1.0; a[1] = 1.0; b[0] = 1.0; b[1] = -1.0;
  double ip(int *, double *, double *); printf("%f \n", ip(&d, a, b));
  double _Complex v[d], w[d]; v[0] = 1; v[1] = 1; w[0] = 1; w[1] = I;
  double _Complex ipc, ip_c(int *, double _Complex *, double _Complex *);
  ipc = ip_c(&d, v, w); printf("%f %f \n", creal(ipc), cimag(ipc));
  return 0;
}
*/


// gcc mat_func.c -lm