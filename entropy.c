//-----------------------------------------------------------------------------------------------------------------------------------
//#include <stdio.h>
//#include <math.h>
//-----------------------------------------------------------------------------------------------------------------------------------
double neumann(int* d, double _Complex rho[][*d]) //! Returns the von Neumann entropy of a density matrix
{ // d  ! Dimension of the density matrix
  // rho(d,d)  ! Density matrix
  double _Complex A[*d][*d];
  int j, k;
  for(j = 0; j < (*d); j++)
  {
    for(k = 0; k < (*d); k++)
    {
      A[j][k] = rho[j][k];
    }
  }
  char jobz = 'N';
  double W[*d];
  int lapacke_zheevd();  lapacke_zheevd(&jobz, d, A, W);
  double vNE, shannon();  // Variables for the shannon and von Neumann entropies
  vNE = shannon(d, W);
  return vNE;
}
//-----------------------------------------------------------------------------------------------------------------------------------
double shannon(int* d, double* pv) // Returns the Shannon entropy of a probability vector
{ // d  ! Dimension of the probability vector
  // pv(d)  ! probability vector
  int j;  // Auxiliary variable for counters
  double SE = 0.0;
  double log2();
  for(j = 0; j < (*d); j++)
  { 
    if( (pv[j] > 1.e-15) && (pv[j] < (1.0-1.e-15)) ) 
    { 
      SE += pv[j]*log2(pv[j]);
    }  
  }
  if (SE > 0.0) SE = -SE;
  return SE;
}
//-----------------------------------------------------------------------------------------------------------------------------------