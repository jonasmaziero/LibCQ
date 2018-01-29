//-----------------------------------------------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
//-----------------------------------------------------------------------------------------------------------------------------------
void coh_test(){
  double theta = 0.0;
  double phi = 0.0;
  double delta = 0.1;
  double pi = M_PI;
  int d; d = 2;
  double _Complex psi[d];
  double _Complex rho[d][d];
  double coh1, coh2, coh_l1(), coh_re();
  int psi_1qb();
  theta = 0.0;
  while (theta < pi)
  {
    psi_1qb(&theta, &phi, psi);  // calls the function which modifies psi in memory
    //printf("%f +I*%f, %f +I*%f \n", creal(psi[0]), cimag(psi[0]), creal(psi[1]), cimag(psi[1]));
    int projector();  projector(&d, psi, rho);  // calls the function which modifies rho in memory 
    coh1 = coh_l1(&d, rho);  coh2 = coh_re(&d, rho);  
    printf("%f \t %f \t %f \n", theta, coh1, coh2);
    theta += delta;
  }
}
//-----------------------------------------------------------------------------------------------------------------------------------
double coh_l1(int *d, double _Complex rho[][*d]){  // Returns the l1-norm coherence
// Ref: T. Baumgratz, M. Cramer e M. B. Plenio, Quantifying coherence, PRL 113, 140401 (2014)
  int j, k;
  double coh = 0.0;  
  for(j = 0; j < ((*d)-1); j++){
    for(k = j+1; k < (*d); k++){  // sum only the above-the-diagonal elements
      coh += sqrt( pow(creal(rho[j][k]),2.0) + pow(cimag(rho[j][k]),2.0) ); // pow(double x, double y) = x**y
    }   
  }
  return 2.0*coh;
}
//-----------------------------------------------------------------------------------------------------------------------------------
double coh_re(int *d, double _Complex rho[][*d]){  //Returns the relative entropy of quantum coherence
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