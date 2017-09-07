//-----------------------------------------------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
//-----------------------------------------------------------------------------------------------------------------------------------
int main()
{
  //tLapack();
  //tTrace();
  //tPTrace();
  //tEntropy();
  int tCoherence();  tCoherence();
  return 0;
}
//-----------------------------------------------------------------------------------------------------------------------------------
int tEntropy()
{
  int d = 2;
  //double pv[] = {0.7,0.3};
  //double SE, shannon();  SE = shannon(&d, pv);  printf("%f \n", SE); // ok
  double _Complex rho[d][d]; rho[0][0] = 1.0;  rho[0][1] = 0.0;  rho[1][0] = 0.0;  rho[1][1] = 0.0;
  double vNE, neumann();  vNE = neumann(&d, rho);  printf("%f \n", vNE);
  return 0;
}
//-----------------------------------------------------------------------------------------------------------------------------------
int tLapack()
{
  int D = 2;
  char jobz = 'V';
  double _Complex B[D][D];
  B[0][0] = 1.0/sqrt(2.0); B[0][1] = 1.0/sqrt(2.0); B[1][0] = 1.0/sqrt(2.0); B[1][1] = -1.0/sqrt(2.0);
  double Ev[D];
  int lapacke_zheevd();  lapacke_zheevd(&jobz, &D, B, Ev);
  printf("%f,%f \n \n", Ev[0], Ev[1]); 
  printf("%f + %f*I, \t %f + %f*I \n", creal(B[0][0]), cimag(B[0][0]), creal(B[1][0]), cimag(B[1][0]));
  printf("%f + %f*I, \t %f + %f*I \n", creal(B[0][1]), cimag(B[0][1]), creal(B[1][1]), cimag(B[1][1]));
  return 0;
}
//-----------------------------------------------------------------------------------------------------------------------------------
int tCoherence()
{
  double theta = 0.0;
  double phi = 0.0;
  double delta = 0.2;
  int d; d = 2;
  double _Complex psi[d];
  double _Complex rho[d][d];
  double coh1, coh2, coh_l1(), coh_re();
  int psi_1qb();
  //while(theta < 2.0*M_PI)
  //{
  theta = M_PI;
    psi_1qb(&theta, &phi, psi);  // calls the function which modifies psi in memory
    printf("%f +I*%f, %f +I*%f \n", creal(psi[0]), cimag(psi[0]), creal(psi[1]), cimag(psi[1]));
    int projector();  projector(&d, psi, rho);  // calls the function which modifies rho in memory
    printf("%f +I*%f, %f +I*%f \n", creal(rho[0][0]), cimag(psi[0]), creal(psi[1]), cimag(psi[1]));
    
    coh1 = coh_l1(&d, rho);
    coh2 = coh_re(&d, rho);  
    printf("%f \t %f \t %f \n", theta, coh1, coh2);
    theta += delta;    
  //}
  return 0;
}
//-----------------------------------------------------------------------------------------------------------------------------------
int tTrace(){
  int N = 2;
  double _Complex HM[N][N];
  double trace_he(), tr_he;
  double _Complex trace_ge(), tr_ge;
  HM[0][0] = 1.0;
  HM[0][1] = 1.0;
  HM[1][0] = 1.0;
  HM[1][1] = 1.0;
  tr_he = trace_he(&N,HM);
  printf("%10.5f \n", tr_he);
  tr_ge = trace_ge(&N,HM);
  printf("%10.5f +I*%10.5f \n", creal(tr_ge), cimag(tr_ge));
  return 0;}
//-----------------------------------------------------------------------------------------------------------------------------------
int tPTrace(){
  int da = 2, db = 2, d = da*db;
  double _Complex rho[d][d];  // In C, arrays are automatically initialized to zero
  int j, k;
  double _Complex rho_b[db][db];
  
  rho[0][0] = 1.0/2.0;  // phi+ Bell state
  rho[0][3] = 1.0/2.0;
  rho[3][0] = 1.0/2.0;
  rho[3][3] = 1.0/2.0;
  int partial_trace_a();  partial_trace_a(&da, &db, rho, rho_b);
  for (j = 0; j < db ; j++)
  {
    for (k = 0; k < db ; k++)
    {
      printf("%f + %f*I \n", creal(rho_b[j][k]), cimag(rho_b[j][k]));
    }
  }
  return 0;}
//-----------------------------------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------------------------------