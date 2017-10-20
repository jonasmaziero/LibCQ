//-----------------------------------------------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
//-----------------------------------------------------------------------------------------------------------------------------------
void main()
{
  //tLapack();
  //tTrace();
  void test_partial_trace(); test_partial_trace();
  //tEntropy();
  //int tCoherence();  tCoherence();
  //tKroneckerP();
  //void test_rng();  test_rng();
  //return 0;
}
//-----------------------------------------------------------------------------------------------------------------------------------
/*void test_rng()
{
  int j;
  FILE *fp;
  char output[] = "output.txt";
  unsigned long long int seed;
  double genrand64_real1();
  srand(time(NULL));  seed = rand();  init_genrand64(seed);
  fp = fopen(output, "w+");
  for (j = 0; j < 10001; j++)
  {
    fprintf(fp, "%15.10f \t %15.10f \n", genrand64_real1(), genrand64_real1());
  }
  for (j = 0; j < 10; j++)
  {
    fscanf(fp, "%15.10f \t %15.10f \n", genrand64_real1(), genrand64_real1());
    printf("%15.10f \t %15.10f \n", genrand64_real1(), genrand64_real1());
  }
  fclose(fp);
}
//-----------------------------------------------------------------------------------------------------------------------------------
int tKroneckerP()
{
  int ra, ca, rb, cb, r, c; 
  ra = 2; ca = 2; rb = 2; cb = 2; r = ra*rb; c = ca*cb;
  double _Complex A[ra][ca], B[rb][cb], KP[r][c];
  A[0][0] = 1.0; A[0][1] = 2.0; A[1][0] = 3.0; A[1][1] = 4.0;
  B[0][0] = 1.0; B[0][1] = 2.0; B[1][0] = 3.0; B[1][1] = 4.0;
  int kronecker_product(), j, k;  
  j = kronecker_product(&ra, &ca, &rb, &cb, A, B, KP);
  for (k = 0; k < r; k++)
  {
    printf("%f %f %f %f \n", creal(KP[k][0]), creal(KP[k][1]), creal(KP[k][2]), creal(KP[k][3])); 
  }
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
  double trace(), tr;
  double _Complex trace_ge(), tr_ge;
  HM[0][0] = 1.0;
  HM[0][1] = 1.0;
  HM[1][0] = 1.0;
  HM[1][1] = 1.0;
  tr = trace(&N,HM);
  printf("%10.5f \n", tr);
  tr_ge = trace_ge(&N,HM);
  printf("%10.5f +I*%10.5f \n", creal(tr_ge), cimag(tr_ge));
  return 0;
}*/
//-----------------------------------------------------------------------------------------------------------------------------------
void test_partial_trace()
{
  int da = 2, db = 4, d = da*db;
  double _Complex rho[d][d];  // In C, arrays are automatically initialized to zero
  int j, k;
  double _Complex rho_b[db][db], rho_a[da][da];
  
  // phi+ Bell state
  //rho[0][0] = 1.0/2.0;  rho[0][3] = 1.0/2.0;  rho[3][0] = 1.0/2.0;  rho[3][3] = 1.0/2.0;
  // W state
  rho[1][1] = 1.0/3.0;  rho[1][2] = 1.0/3.0;  rho[1][4] = 1.0/3.0;  
  rho[2][1] = 1.0/3.0;  rho[2][2] = 1.0/3.0;  rho[2][4] = 1.0/3.0; rho[4][1] = 1.0/3.0;  rho[4][2] = 1.0/3.0;  rho[4][4] = 1.0/3.0;
  partial_trace_a(&da, &db, rho, rho_b);
  //partial_trace_b(&da, &db, rho, rho_a);
  array_display(&db, &db, rho_b);
}
//-----------------------------------------------------------------------------------------------------------------------------------