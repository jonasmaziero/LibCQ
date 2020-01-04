#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>

double ip(int *d, double *v, double *w) {
  int j;
  double inner = 0.0;
  for(j = 0; j < (*d); j++) {
    if (fabs(*(v+j)) > 1.e-15 && fabs(*(w+j)) > 1.e-15) { 
      inner += (*(v+j))*(*(w+j));
    }
  }
  return inner;
}

double norm(int *d, double *v) {
  double ip(int *, double *, double *);  
  return sqrt(ip(d, v, v));
}

double _Complex ip_c(int *d, double _Complex *v, double _Complex *w) {
  int j;
  double _Complex inner = 0.0;
  for (j = 0; j < (*d); j++) {
    inner += creal(*(v+j))*creal(*(w+j)) + cimag(*(v+j))*cimag(*(w+j));
    inner += I*(creal(*(v+j))*cimag(*(w+j)) - cimag(*(v+j))*creal(*(w+j)));
  }
  return inner;
}

double norm_c(int *d, double _Complex *v){
  double nm = 0.0;
  double _Complex ip_c(int *, double _Complex *, double _Complex *);
  double _Complex inner = ip_c(d, v, v);
  nm = sqrt(pow(creal(inner),2.0) + pow(cimag(inner),2.0));
  return nm;
}

/*
int main() {
  int d = 2;
  double *v; v = (double *)malloc(d*sizeof(double));
  double _Complex *w; w = (double _Complex *)malloc(d*sizeof(double _Complex));
  *(v+0) = -1.0/sqrt(2.0);  *(v+1) = 1.0/sqrt(2.0);  
  *(w+0) = 1.0/sqrt(2.0);  *(w+1) = -I/sqrt(2.0);
  double norm(int *, double *);  printf("%f \n", norm(&d, v));
  double norm_c(int *, double _Complex *);  printf("%f \n", norm_c(&d, w));
  free(v); free(w);
  return 0;
}
*/