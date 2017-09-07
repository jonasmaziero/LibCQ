//-----------------------------------------------------------------------------------------------------------------------------------
//#include <stdio.h>
#include <complex.h>
//-----------------------------------------------------------------------------------------------------------------------------------
double trace_he(int* N, double _Complex HM[][*N]) // Computes the trace of an HERMITIAN NxN matrix
{
  int j;
  double trace = 0.0; 
  for(j = 0; j < (*N); j++)
  { 
    trace += creal(HM[j][j]);
  }    
  return trace;
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
int projector(int* d, double _Complex* psi, double _Complex rho[][*d])
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
  return 0;
}
//-----------------------------------------------------------------------------------------------------------------------------------