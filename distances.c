//-----------------------------------------------------------------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <complex.h>
//-----------------------------------------------------------------------------------------------------------------------------------
void innerTest(void){
  int d = 2;
  //double v[d], w[d];
  //double _Complex v[d], w[d];
  double v[d];
  double _Complex w[d];
  v[0] = -1.0/sqrt(2.0);  v[1] = 1.0/sqrt(2.0);  w[0] = 1.0/sqrt(2.0);  w[1] = -I/sqrt(2.0);
  //double innerR();  printf("%f \n", innerR(&d, v, w)));
  //double _Complex innerC();  printf("%f %f \n", creal(innerC(&d, v, w)), cimag(innerC(&d, v, w)));
  double normR();  printf("%f \n", normR(&d, v));
  double normC();  printf("%f \n", normC(&d, w));
}
//-----------------------------------------------------------------------------------------------------------------------------------
double innerR(int *d, double *v, double *w){
  int j;
  double ip = 0.0;
  for(j = 0; j < (*d); j++){
    if (fabs(v[j]) > 1.e-15 && fabs(w[j]) > 1.e-15) ip += (v[j])*(w[j]);
  }
  return ip;
}
//-----------------------------------------------------------------------------------------------------------------------------------
double _Complex innerC(int *d, double _Complex *v, double _Complex *w){
  int j;
  double _Complex ip = 0.0;
  for(j = 0; j < (*d); j++){
    ip = ip + creal(v[j])*creal(w[j]) + cimag(v[j])*cimag(w[j]) + I*(creal(v[j])*cimag(w[j]) - cimag(v[j])*creal(w[j]));
  }
  return ip;
}
//-----------------------------------------------------------------------------------------------------------------------------------
double normR(int *d, double *v){
  double norm = 0.0;
  double innerR();  norm = innerR(d, v, v);
  norm = sqrt(norm);
  return norm;
}
//-----------------------------------------------------------------------------------------------------------------------------------
double normC(int *d, double _Complex *v){
  double norm = 0.0;
  double _Complex innerC(), ip;  ip = innerC(d, v, v);
  norm = sqrt(pow(creal(ip),2.0) + pow(cimag(ip),2.0));
  return norm;
}
//-----------------------------------------------------------------------------------------------------------------------------------
