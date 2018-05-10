//------------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <complex.h>
//------------------------------------------------------------------------------
void idR(int *d, double id[][*d]){
  // returns the identity matrix of dimension d
  int j, k;
  for(j = 0; j < (*d); j++){
    for(k = j; k < (*d); k++){
      if(j = k){id[j][k] = 1.0;} else {id[j][k] = 0.0;  id[k][j] = 0.0;}
    }
  }
}
//------------------------------------------------------------------------------
void idC(int *d, double _Complex id[][*d]){
  // returns the identity matrix of dimension d
  int j, k;
  for(j = 0; j < (*d); j++){
    for(k = j; k < (*d); k++){
      if(j = k){id[j][k] = 1.0;} else {id[j][k] = 0.0;  id[k][j] = 0.0;}
    }
  }
}
//------------------------------------------------------------------------------
double maxArray1D(int *d, double *A){
  // returns the maximum value among those stored in a double precision 1D array
  int j;  double max = A[0];
  for(j = 1; j < (*d); j++){if(A[j] > max)max = A[j];}
  return max;
}
//------------------------------------------------------------------------------
double traceC(int *d, double _Complex *A){
  // Computes the trace of a general complex NxN matrix
  int j;
  double tr = 0.0;
  for(j = 0; j < (*d); j++){
    tr += (*(A+j*(*d)+j));
  }
  return tr;
}
//------------------------------------------------------------------------------
double traceHe(int *N, double _Complex *HM){
  // Computes the trace of an HERMITIAN NxN matrix
  int j;
  double tr = 0.0;
  for(j = 0; j < (*N); j++){
    tr += creal(*(HM+j*(*N)+j));
  }
  return tr;
}
//------------------------------------------------------------------------------
double _Complex trace_ge(int *N, double _Complex GM[][*N]){
  // Computes the trace of a GENERAL NxN matrix
  int j;
  double _Complex trace = 0.0;
  for(j = 0; j < (*N); ++j){
    trace += GM[j][j];
  }
  return trace;
}
//------------------------------------------------------------------------------
void projector(int *d, double _Complex *psi, double _Complex rho[][*d]){
  // returns the projector on state psi
  // d  ! Dimension of the vector
  // psi(d)  ! Vector we want the projector on
  // rho(d,d)  ! Projector on vec
  int j, k;
  for(j = 0; j < (*d); j++){
    for(k = j; k < (*d); k++){
      rho[j][k] = (psi[j])*(conj(psi[k]));
      // Elements in the diagonal and above
      if(j != k){ rho[k][j] = conj(rho[j][k]); }
      // Elements below the diagonal
    }
  }
}
//------------------------------------------------------------------------------
void kronecker(int *ra, int *ca, int *rb, int *cb,
                double _Complex A[][*ca], double _Complex B[][*cb],
                double _Complex kp[][(*ca)*(*cb)]){
  // Returns the Kronecker product of two general complex matrices
  int ja, ka, jb, kb, j, k, jl, ju, kl, ku;
  for (ja = 0; ja < (*ra); ja++){
    jl = ja*(*rb);  ju = jl+(*rb);
    for (ka = 0; ka < (*ca); ka++){
      if (A[ja][ka] != 0.0){
        jb = -1;
        kl = ka*(*cb);  ku = kl+(*cb);
        for (j = jl; j < ju; j++){
          jb += 1;
          kb = -1;
          for (k = kl; k < ku; k++){
            kb += 1;
            kp[j][k] = A[ja][ka]*B[jb][kb];
            /*kp[j][k] = creal(A[ja][ka])*creal(B[jb][kb])
                         - cimag(A[ja][ka])*cimag(B[jb][kb])
                         + I*(creal(A[ja][ka])*cimag(B[jb][kb])
                         + cimag(A[ja][ka])*creal(B[jb][kb])); */
          }
        }
      }
    }
  }
}
//------------------------------------------------------------------------------
void matmulR(int *ra, int *ca, int *cb, double A[][*ca], double B[][*cb],
             double C[][*cb]) {
  // returns the product of two matrices
  int j, k, l;
  for(j = 0; j < (*ra); j++){
    for(k = 0; k < (*cb); k++){
      C[j][k] = 0.0;
      for(l = 0; l < (*ca); l++){
        C[j][k] += A[j][l]*B[l][k];
      }
    }
  }
}
//------------------------------------------------------------------------------
void matmulC(int *ra, int *ca, int *cb, double _Complex A[][*ca],
             double _Complex B[][*cb], double _Complex C[][*cb]) {
  // returns the product of two matrices
  int j, k, l;
  for(j = 0; j < (*ra); j++){
    for(k = 0; k < (*cb); k++){
      C[j][k] = 0.0;
      for(l = 0; l < (*ca); l++){
        C[j][k] += creal(A[j][l])*creal(B[l][k])
                   - cimag(A[j][l])*cimag(B[l][k])
                   + I*(creal(A[j][l])*cimag(B[l][k])
                   + cimag(A[j][l])*creal(B[l][k]));
      }
    }
  }
}
//------------------------------------------------------------------------------
void array2DisplayR(int *xd, int *yd, double *A){
  // Shows an 2D real array in the screen
  int j, k;
  for(j = 0; j < (*xd); j++){
    for(k = 0; k < (*yd); k++){
      printf("%f \t", *(A+k+j*(*yd)));
    }
    printf("\n");
  }
}
//------------------------------------------------------------------------------
void array2DisplayC(int *xd, int *yd, double _Complex *A){
  // Shows the real and imaginary parts of an array in the screen
  int j, k;
  printf("real part \n");
  for (j = 0; j < (*xd); j++){
    for (k = 0; k < (*yd); k++){
      printf("%lf \t", creal(*(A+k+j*(*yd))));
    }
    printf("\n");
  }
  printf("imaginary part \n");
  for (j = 0; j < (*xd); j++){
    for (k = 0; k < (*yd); k++){
      printf("%lf \t", cimag(*(A+k+j*(*yd))));
    }
    printf("\n");
  }
}
//------------------------------------------------------------------------------
void transposeR(int *nr, int *nc, double A[][*nc], double At[][*nr]){
  int j, k;
  for(j = 0; j < (*nr); j++){
    for(k = 0; k < (*nc); k++){
      At[k][j] = A[j][k];
    }
  }
}
//------------------------------------------------------------------------------
void transposeC(int *nr, int *nc, double _Complex A[][*nc], double _Complex At[][*nr]){
  int j, k;
  for(j = 0; j < (*nr); j++){
    for(k = 0; k < (*nc); k++){
      At[k][j] = A[j][k];
    }
  }
}
//------------------------------------------------------------------------------
void adjoint(int *nr, int *nc, double _Complex A[][*nc], double _Complex Ad[][*nr]){
  int j, k;
  for(j = 0; j < (*nr); j++){
    for(k = 0; k < (*nc); k++){
      Ad[k][j] = creal(A[j][k]) - cimag(A[j][k])*I;
    }
  }
}
//------------------------------------------------------------------------------
void matconj(int *nr, int *nc, double _Complex A[][*nc], double _Complex Ac[][*nr]){
  int j, k;
  for(j = 0; j < (*nr); j++){
    for(k = 0; k < (*nc); k++){
      Ac[j][k] = creal(A[j][k]) - cimag(A[j][k])*I;
    }
  }
}
//------------------------------------------------------------------------------
void matA2Bc(int *xd, int *yd, double _Complex *A, double _Complex *B) {
  int j, k;
  for (j = 0; j < (*xd); j++) {
    for (k = 0; k < (*yd); k++) {
      *(B+j*(*yd)+k) = *(B+j*(*yd)+k);
    }
  }
}
//------------------------------------------------------------------------------
