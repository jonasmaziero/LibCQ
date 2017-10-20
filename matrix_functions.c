//-----------------------------------------------------------------------------------------------------------------------------------
#include <stdio.h>
#include <complex.h>
//-----------------------------------------------------------------------------------------------------------------------------------
double trace(int* N, double _Complex HM[][*N]) // Computes the trace of an HERMITIAN NxN matrix
{
  int j;
  double tr = 0.0; 
  for(j = 0; j < (*N); j++)
  { 
    tr += creal(HM[j][j]);
  }    
  return tr;
}
//-----------------------------------------------------------------------------------------------------------------------------------
double _Complex trace_ge(int* N, double _Complex GM[][*N]) // Computes the trace of a GENERAL NxN matrix
{
  int j;
  double _Complex trace = 0.0; 
  for(j = 0; j < (*N); ++j)
  { 
    trace += GM[j][j];
  }  
  return trace;
}
//-----------------------------------------------------------------------------------------------------------------------------------
void projector(int* d, double _Complex* psi, double _Complex rho[][*d])  // returns the projector on state psi
{ // d  ! Dimension of the vector
  // psi(d)  ! Vector we want the projector on
  // rho(d,d)  ! Projector on vec
  int j, k;
  for(j = 0; j < (*d); j++)
  {
    for(k = j; k < (*d); k++)
    {
      rho[j][k] = (psi[j])*(conj(psi[k])); // Elements in the diagonal and above
      if(j != k){ rho[k][j] = conj(rho[j][k]); } // Elements below the diagonal
    }  
  }
}
//------------------------------------------------------------------------------------------------------------------------------------
// Returns the tensor product of two general complex matrices
void kronecker_product(int* ra, int* ca, int* rb, int* cb, 
                      double _Complex A[][*ca], double _Complex B[][*cb], double _Complex KP[][(*ca)*(*cb)])  
{
  // ra, ca, rb, cb  ! Number of rows and columns of the two matrices
  // M1(ra,ca), M2(rb,cb)  ! Matrices to take the tensor product of
  // M1_kp_M2(ra*rb,ca*cb)  ! Matrix containing the tensor product of M1 and M2
  int ja, ka, jb, kb, j, k, jl, ju, kl, ku;  // Auxiliary variables for counters
  for (ja = 0; ja < (*ra); ja++)
  {
    jl = ja*(*rb);  ju = jl+(*rb);
    for (ka = 0; ka < (*ca); ka++)
    {
      if (A[ja][ka] != 0.0)
      {
        jb = -1;
        kl = ka*(*cb);  ku = kl+(*cb);
        for (j = jl; j < ju; j++)
        { 
          jb += 1;
          kb = -1;
          for (k = kl; k < ku; k++)
          {
            kb += 1;
            KP[j][k] = A[ja][ka]*B[jb][kb];
          }
        }
      }
    }
  }
}
//-----------------------------------------------------------------------------------------------------------------------------------
// Shows the real and imaginary parts of an array in the screen
void array_display(int* r, int* c, double _Complex array[][*c])
{
  int j, k;
  printf("real part \n");
  for (j = 0; j < (*r); j++)
  {
    for (k = 0; k < (*c); k++)
    {
      printf("%f \t", creal(array[j][k]));
    }
    printf("\n");
  }
  printf("imaginary part \n");
  for (j = 0; j < (*r); j++)
  {
    for (k = 0; k < (*c); k++)
    {
      printf("%f \t", cimag(array[j][k]));
    }
    printf("\n");
  }
}
//-----------------------------------------------------------------------------------------------------------------------------------