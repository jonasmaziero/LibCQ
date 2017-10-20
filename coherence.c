//-----------------------------------------------------------------------------------------------------------------------------------
#include <stdlib.h>
#include <math.h>
#include <complex.h>
//-----------------------------------------------------------------------------------------------------------------------------------
double coh_l1(int* d, double _Complex rho[][*d]){  // Returns the l1-norm coherence
// Ref: T. Baumgratz, M. Cramer e M. B. Plenio, Quantifying coherence, PRL 113, 140401 (2014)
  int j, k;
  double coh = 0.0;  
  for(j = 0; j < ((*d)-1); j++){
    for(k = j+1; k < (*d); k++){  // sum only the above-the-diagonal elements
      coh += sqrt( pow(2,creal(rho[j][k])) + pow(2,cimag(rho[j][k])) );
    }   
  } 
  coh = 2.0*coh;
  return coh;
}
//-----------------------------------------------------------------------------------------------------------------------------------
double coh_re(int* d, double _Complex rho[][*d]){  //Returns the relative entropy of quantum coherence
// Ref: T. Baumgratz, M. Cramer e M. B. Plenio, Quantifying coherence, PRL 113, 140401 (2014)
  int j;
  double pv[*d];
  for(j = 0; j < (*d); j++){
    pv[j] = creal(rho[j][j]);   
  }
  double coh, shannon(), neumann();  
  coh = shannon(d, pv) - neumann(d, rho);
  return coh;  
}
//-----------------------------------------------------------------------------------------------------------------------------------