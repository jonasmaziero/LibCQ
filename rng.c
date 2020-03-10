#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

void rng_init() { 
  unsigned long long int seed;  void srand();  srand(time(NULL));
  int rand();  seed = rand(); void init_genrand64();  init_genrand64(seed);
}

void rng_gauss(double *grn1, double *grn2) {
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

void rng_gaussl(long double *grn1, long double *grn2) {
  long double u, v, s = 2.0, w, a = -1.0, b = 1.0;
  void rn_abl(long double *, long double *, long double *);
  do {
    rn_abl(&a, &b, &u); rn_abl(&a, &b, &v); s = powl(u,2)+powl(v,2);
  } while ( s >= 1 || s == 0);
  w = sqrtl(-2.0*logl(s)/s);
  (*grn1) = u*w; (*grn2) = v*w;
}

void rng_exp(double *ern) {
  double genrand64_real1();
  double rn = 0.0;
  while (rn == 0.0) rn = genrand64_real1();
  *(ern) = -log(rn);
}

void rng_expl(long double *ern) { // long double version
  long double genrand64_real1l();
  long double rn = 0.0;
  while (rn == 0.0) rn = genrand64_real1l();
  *(ern) = -logl(rn);
}

void rn_ab(double *a, double *b, double *rn) {
  double genrand64_real1();
  (*rn) = (*a) + ((*b)-(*a))*genrand64_real1();
}

void rn_abl(long double *a, long double *b, long double *rn) {
  long double genrand64_real1l();
  (*rn) = (*a) + ((*b)-(*a))*genrand64_real1l();
}

int rn_int_ab(int *a, int *b) {
  long long genrand64_int63();
  int ri;
  ri = ((*a) + ((*b)-(*a)+1)*(genrand64_int63()/(pow(2,63)-1)));
  return ri;
}

/*
void test_gauss() {
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
  for (l = 0; l < (ni-1); l++) { s += ct[l]; }; printf("%f \n",(double)s/(2*(double)ns));
  FILE *fd = fopen("plot.dat", "w");
  for(l = 0; l < (ni-1); l++){
    fprintf(fd, "%f %f \n", (xmin +(l+0.5)*delta), (double)ct[l]/(delta*(double)(2*ns)));
  }
  fclose(fd);
  void plot(); plot();
}

void test_exp() {
  int j, l, ns = pow(10,4), ni = 200;
  double xmin = 0.0, xmax = 5.0;
  double delta = (xmax-xmin)/((double) ni);
  int ct[ni];
  for(j = 0; j < ni; j++){
    ct[j] = 0;
  }
  void rng_exp(double *); 
  double ern;
  for(j = 0; j < ns; j++){
    rng_exp(&ern);
    for (l = 0; l < (ni-1); l++) {
      if(ern >= (xmin + l*delta) && ern < (xmin + (l+1)*delta)) {
        ct[l] += 1;
      }
    }
  }
  int s = 0;
  for (l = 0; l < (ni-1); l++) { s += ct[l]; }; 
  printf("%f \n",((double) s)/((double) ns));
  FILE *fd = fopen("plot.dat", "w");
  for(l = 0; l < (ni-1); l++){
    fprintf(fd, "%f %f \n", (xmin +(l+0.5)*delta), ((double) ct[l])/(delta*ns));
  }
  fclose(fd);
  void plot(); plot();
}

void test_rand_int() {
  int j, ns = pow(10,4), a, b, ri, rn_int_ab(int *, int *);
  int *ct; ct = (int *)malloc((b-a+1)*sizeof(int));
  double pd;
  a = 3; b = 10;
  for (j = 0; j < ns; j++) {
    ri = rn_int_ab(&a, &b);
    *(ct+ri) += 1;
  }
  FILE *fd = fopen("plot.dat", "w");
  for (j = 0; j < (b-a+1); j++) {
    pd = ((double) *(ct+j+a))/((double) ns);
    fprintf(fd, "%i %f \n", j+a, pd);
  }
  fclose(fd);
  free(ct);
  void plot(); plot();
}

void six_int() {
  int rn_int_ab(int *, int *); 
  int a, b, j, k, ni = 0, ri, si[6], re; 
  a = 0; b = 60;
  ri = rn_int_ab(&a, &b); si[0] = ri; ni = 0;
  while (ni < 6) {
    ri = rn_int_ab(&a, &b);
    re = 0;
    for (k = 0; k < (ni+1); k++) {
      if (ri == si[k]) re = 1;
    }
    if (re == 0) { ni += 1; si[ni] = ri; }
  }
  for (j = 0; j < 6; j++) { printf("%i \t", si[j]); }
  printf("\n");
}


int main() {
  //void rng_init(); rng_init();
  //void test_gauss(); test_gauss();
  //void test_exp(); test_exp();
  //void test_rand_int(); test_rand_int();
  //void six_int(); six_int();
  long double genrand64_real1l(); printf("%.20Lf \n", genrand64_real1l());
  return 0;
}
*/
// gcc rng.c MT19937_64.c plots.c -lm

