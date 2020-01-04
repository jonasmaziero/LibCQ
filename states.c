
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

void psi1qb(double *theta, double *phi, double _Complex *psi) {
  *(psi+0) = cos((*theta)/2.0);  
  *(psi+1) = sin((*theta)/2.0)*(cos((*phi)) + I*sin((*phi)));
}

void bell(double _Complex *phip, double _Complex *psip, double _Complex *phim, double _Complex *psim){
  double f = 1.0/sqrt(2.0);
  phip[0] = f; phip[1] = 0.0; phip[2] = 0.0; phip[3] = f;
  phim[0] = f; phim[1] = 0.0; phim[2] = 0.0; phim[3] = -f;
  psip[0] = 0.0; psip[1] = f; psip[2] = f; psip[3] = 0.0;
  psim[0] = 0.0; psim[1] = f; psim[2] = -f; psim[3] = 0.0;
}

void werner2qb(double *p, double _Complex rho[][4]){
  double _Complex I4[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
  double _Complex proj[4][4] = {{0,0,0,0},{0,0.5,-0.5,0},{0,-0.5,0.5,0},{0,0,0,0}};
  int j, k;
  for(j = 0; j < 4; j++){
    for(k = 0; k < 4; k++){
      rho[j][k] = (1.0-(*p))*0.25*I4[j][k] + (*p)*proj[j][k];
    }
  }
}
