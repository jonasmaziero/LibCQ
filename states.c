//-----------------------------------------------------------------------------------------------------------------------------------
#include <stdio.h>
//#include <stdlib.h>
#include <math.h>
#include <complex.h>
//-----------------------------------------------------------------------------------------------------------------------------------
int psi_1qb(double* theta, double* phi, double _Complex* psi)
{
  psi[0] = cos((*theta)/2.0);  psi[1] = sin((*theta)/2.0)*(cos((*phi)) + I*sin((*phi)));
  printf("%f +I*%f, %f +I*%f \n", creal(psi[0]), cimag(psi[0]), creal(psi[1]), cimag(psi[1]));
  return 0;
}
//-----------------------------------------------------------------------------------------------------------------------------------