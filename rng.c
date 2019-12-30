// compile with: gcc rng.c MT19937_64.c -lm
//------------------------------------------------------------------------------
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "MT19937_64.h"
//------------------------------------------------------------------------------
void rng_test(){
  unsigned long long int seed;  srand(time(NULL));  seed = rand();
  void init_genrand64();  init_genrand64(seed);
  void test_gauss(void); test_gauss();
}
//------------------------------------------------------------------------------
void test_gauss(){
  printf("Computing the probability distribution\n");
  int j, l;
  int ns = pow(10,4), ni = 100;
  double xmin = -5.0, xmax = 5.0;
  double delta =  (xmax-xmin)/((double) ni);
  int ct[ni];
  for(j = 0; j < ni; j++){
    ct[j] = 0;
  }
  void rng_gauss(double *, double *);
  double grn1, grn2;
  for(j = 0; j < ns; j++){
    rng_gauss(&grn1, &grn2);
    for (l = 0; l < (ni-1); l++) {
      if((grn1 >= (xmin + l*delta)) && (grn1 < (xmin + (l+1)*delta))) {
        ct[l] = ct[l] + 1;
      }
      if((grn2 >= (xmin + l*delta)) && (grn2 < (xmin + (l+1)*delta))) {
        ct[l] = ct[l] + 1;
      }
    }
  }
  int s = 0;
  for (l = 0; l < (ni-1); l++) { s += ct[l]; }; printf("%f\n",(double)s/(2*(double)ns));
  FILE *fd = fopen("rn.dat", "w");
  for(l = 0; l < (ni-1); l++){
    fprintf(fd, "%f %f \n", (xmin +(l+0.5)*delta), (double)ct[l]/(delta*(double)(2*ns)));
  }
  fclose(fd);
  // writing a gnuplot script, from C
  FILE *fgnu = fopen("rn.gnu", "w");
  fprintf(fgnu, "reset \n");
  fprintf(fgnu, "set terminal postscript enhanced 'Helvetica' 24 \n");
  fprintf(fgnu, "set output 'rn.eps' \n");
  fprintf(fgnu, "plot 'rn.dat' t 'numerico',(1/sqrt(2*pi))*exp(-x**2/2)  t 'exato' \n");
  fclose(fgnu);
  system("gnuplot rn.gnu");
  system("evince rn.eps&");
}
//------------------------------------------------------------------------------
void rng_init(){ // initializes the MT random number generator
  unsigned long long int seed;  void srand();  srand(time(NULL));
  int rand();  seed = rand();
  void init_genrand64();  init_genrand64(seed);
}
//------------------------------------------------------------------------------
void rng_gauss(double *grn1, double *grn2){
  // returns random numbers with Gaussian probability distribution
  /*// Box–Muller method
  double genrand64_real1();
  rn = genrand64_real1(); if (rn < 1.e-15) rn = 1.0E-15;
  logterm = sqrt(-2.0*log(rn));  angle = 2.0*M_PI*genrand64_real1();
  (*grn1) = logterm*cos(angle); (*grn2) = logterm*sin(angle); //*/
  // Marsaglia method (avoids sin and cos calcs)
  double u, v, s = 2.0, w, a = -1.0, b = 1.0;
  void rn_ab(double *, double *, double *);
  do {
   rn_ab(&a, &b, &u); rn_ab(&a, &b, &v); s = pow(u,2)+pow(v,2);
 } while ( s >= 1 || s == 0);
  w = sqrt(-2.0*log(s)/s);
  (*grn1) = u*w; (*grn2) = v*w;
}
//------------------------------------------------------------------------------
void rn_ab(double *a, double *b, double *rn){
  // return a random number uniformly distributed in [a,b)
  double genrand64_real1();
  (*rn) = (*a) + ((*b)-(*a))*genrand64_real1();
}
//------------------------------------------------------------------------------
