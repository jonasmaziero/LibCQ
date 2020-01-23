#include <lapacke.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
/* 
install lapacke with: sudo apt-get install liblapacke-dev
To see how to call the subroutines, take a look at the lapacke.h file in 
/usr/local/include
*/

void lapacke_zheevd(char *jobz, int *D, double _Complex B[][*D], double *W) {
  lapack_int info;
  lapack_int d = *D;
  lapack_complex_double A[d][d];
  int j, k;
  for (j = 0; j < (*D); j++) {
    for (k = 0; k < (*D); k++) {
      A[j][k] = B[j][k];
    }
  }
  char JOBZ = *jobz;
  info = LAPACKE_zheevd(LAPACK_ROW_MAJOR, JOBZ, 'U', d, *A, d, W);
  if (*jobz == 'V') {
    for (j = 0; j < (*D); j++) {
      for (k = 0; k < (*D); k++) {
        B[j][k] = A[j][k];
      }
    }
  }
}

void lapacke_zgeev(char *jobz, int *D, double _Complex B[][*D], double _Complex *W) {
  lapack_int info;
  lapack_int d = *D;
  lapack_complex_double Ac[d][d], A[d][d], C[d][d];
  int j, k;
  for (j = 0; j < (*D); j++) {
    for (k = 0; k < (*D); k++) {
      A[j][k] = B[j][k];
    }
  }
  char jobvr = *jobz, jobvl = 'N';
  info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, jobvl, jobvr, d, *Ac, d, W, *C, d, *A, d);
  if (*jobz == 'V') {
    for (j = 0; j < (*D); j++) {
      for (k = 0; k < (*D); k++) {
        B[j][k] = A[j][k];
      }
    }
  }
}

/*
int main() {
  int D = 2;
  char jobz = 'V';
  double _Complex B1[D][D], B2[D][D], egval[D];
  B1[0][0] = 1;  B1[0][1] = 0;  B1[1][0] = 0;  B1[1][1] = -1;
  B2[0][0] = B1[0][0]; B2[0][1] = B1[0][1];
  B2[1][0] = B1[1][0]; B2[1][1] = B1[1][1];
  double Eval[D];
  void lapacke_zheevd(); lapacke_zheevd(&jobz, &D, B1, Eval);
  printf("%f,%f \n", Eval[0], Eval[1]);
  printf("%f + %f*I, \t %f + %f*I \n",
          creal(B1[0][0]), cimag(B1[0][0]), creal(B1[1][0]), cimag(B1[1][0]));
  printf("%f + %f*I, \t %f + %f*I \n",
          creal(B1[0][1]), cimag(B1[0][1]), creal(B1[1][1]), cimag(B1[1][1]));
  void lapacke_zgeev(); lapacke_zgeev(&jobz, &D, B2, egval);
  printf("%f,%f \n", creal(egval[0]), creal(egval[1]));
  printf("%f + %f*I, \t %f + %f*I \n",
          creal(B2[0][0]), cimag(B2[0][0]), creal(B2[1][0]), cimag(B2[1][0]));
  printf("%f + %f*I, \t %f + %f*I \n",
          creal(B2[0][1]), cimag(B2[0][1]), creal(B2[1][1]), cimag(B2[1][1]));
  return 0;
}
*/

//  gcc lapack.c -llapacke