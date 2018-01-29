//-----------------------------------------------------------------------------------------------------------------------------------
#include <stdio.h>
//#include <math.h>
//-----------------------------------------------------------------------------------------------------------------------------------
int entropy_test(){
  int d = 2;
  //double pv[] = {0.7,0.3};
  //double SE, shannon();  SE = shannon(&d, pv);  printf("%f \n", SE); // ok
  double _Complex rho[d][d]; rho[0][0] = 0.5;  rho[0][1] = 0.0;  rho[1][0] = 0.0;  rho[1][1] = 0.5;
  double vNE, neumann();  vNE = neumann(&d, rho);  printf("%f \n", vNE);
  return 0;
}
//-----------------------------------------------------------------------------------------------------------------------------------
double neumann(int *d, double _Complex rho[][*d]){ //! Returns the von Neumann entropy of a density matrix
  // d  ! Dimension of the density matrix
  // rho(d,d)  ! Density matrix
  double _Complex A[*d][*d];
  int j, k;
  for(j = 0; j < (*d); j++){
    for(k = 0; k < (*d); k++){
      A[j][k] = rho[j][k];
    }
  }
  char jobz = 'N';  double W[*d];  lapacke_zheevd(&jobz, d, A, W);
  double vNE, shannon();  // Variables for the shannon and von Neumann entropies
  vNE = shannon(d, W);
  return vNE;
}
//-----------------------------------------------------------------------------------------------------------------------------------
double shannon(int *d, double *pv){ // Returns the Shannon entropy of a probability vector
  // d  ! Dimension of the probability vector
  // pv(d)  ! probability vector
  int j;
  double SE = 0.0;
  double log2();
  for(j = 0; j < (*d); j++){ 
    if((pv[j] > 1.e-15) && (pv[j] < (1.0-1.e-15))){ 
      SE -= pv[j]*log2(pv[j]);
    }  
  }
  return SE;
}
//-----------------------------------------------------------------------------------------------------------------------------------