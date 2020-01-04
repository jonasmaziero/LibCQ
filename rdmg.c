#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>

void rdm_ginibre(int *d, double _Complex *rdm) {
  int j,k,l;
  double _Complex *G; 
  G = (double _Complex *)malloc((*d)*(*d)*sizeof(double _Complex));
  void ginibre(int *, double _Complex *); ginibre(d, G);
  double N2, norm_hs(int *, double _Complex *);
  //void zm_c(int *d, double _Complex *z); zm_c(d, rdm);
  N2 = pow(norm_hs(d,G),2.0);
  for (j = 0; j < (*d); j++) {
    for (k = j; k < (*d); k++) {
      for (l = 0; l < (*d); l++) {
        *(rdm+j*(*d)+k) += (creal(*(G+j*(*d)+l)))*(creal(*(G+k*(*d)+l)));
        *(rdm+j*(*d)+k) += (cimag(*(G+j*(*d)+l)))*(cimag(*(G+k*(*d)+l)));
        *(rdm+j*(*d)+k) -= I*(creal(*(G+j*(*d)+l)))*(cimag(*(G+k*(*d)+l)));
        *(rdm+j*(*d)+k) += I*(cimag(*(G+j*(*d)+l)))*(creal(*(G+k*(*d)+l)));
      }
      *(rdm+j*(*d)+k) /= N2;
      if (j != k) {
        *(rdm+k*(*d)+j) = creal(*(rdm+j*(*d)+k)) - I*cimag(*(rdm+j*(*d)+k));
      }
    }
  }
}

void ginibre(int *d, double _Complex *G){
  double grn1, grn2;
  void rng_gauss(double *, double *);
  int j, k;
  for (j = 0; j < (*d); j++) {
    for (k = 0; k < (*d); k++) {
      rng_gauss(&grn1, &grn2);
      *(G+j*(*d)+k) = grn1 + I*grn2;
    }
  }
}

void rdm_pos_sub(int *d, double _Complex *rrho) {
  double *rpv; 
  rpv = (double *)malloc((*d)*sizeof(double));
  //void rpv_zhsl(int *, double *); rpv_zhsl(d, rpv);
  void rpv_devroye(int *, double *); rpv_devroye(d, rpv);
  void rand_circle(long double *, long double *, long double *);
  long double r, x, y;
  int j, k;
  for (j = 0; j < (*d); j++) {
    *(rrho+j*(*d)+j) = *(rpv+j);
  }
  free(rpv);
  for (j = 0; j < ((*d)-1); j++) {
    for (k = (j+1); k < (*d); k++) {
      r = sqrt(cabs(*(rrho+j*(*d)+j))*cabs(*(rrho+k*(*d)+k)));
      rand_circle(&r, &x, &y);
      *(rrho+j*(*d)+k) = x + I*y; *(rrho+k*(*d)+j) = x - I*y;
    }
  }
}

void rdm_pos(int *d, double _Complex *rrho) {
  double *rpv; 
  rpv = (double *)malloc((*d)*sizeof(double));
  //void rpv_zhsl(int *, double *); rpv_zhsl(d, rpv);
  void rpv_devroye(int *, double *); rpv_devroye(d, rpv);
  void rand_circle(long double *, long double *, long double *);
  long double r, x, y;
  int j, k;
  for (j = 0; j < (*d); j++) {
    *(rrho+j*(*d)+j) = *(rpv+j);
  }
  free(rpv);
  r = 0.0;
  for (j = 0; j < ((*d)-1); j++) {
    for (k = (j+1); k < (*d); k++) {
      r += (*(rrho+j*(*d)+j))*(*(rrho+k*(*d)+k));
    }
  }
  r = sqrt(r);
  for (j = 0; j < ((*d)-1); j++) {
    for (k = j+1; k < (*d); k++) {
      rand_circle(&r, &x, &y);
      *(rrho+j*(*d)+k) = x + I*y;
      *(rrho+k*(*d)+j) = x - I*y;
      r -= (pow(x,2)+pow(y,2));
    }
  }
}

void rand_circle(long double *r, long double *x, long double *y) {
  // returns a random point (x,y) uniformly distributed
  // within a circle of radius r
  long double th, ph, rh;
  double genrand64_real1();
  th = 2.0*M_PI*genrand64_real1();
  rh = (*r)*sqrt(genrand64_real1());
  (*x) = rh*cos(th); 
  (*y) = rh*sin(th);
}

void test_rand_circle() {  // there is a problem for r < 0.00001
  int n, j;
  long double r, x, y;
  void rand_circle(long double *, long double *, long double *);
  FILE *fd = fopen("plot.dat", "w");
  n = 10000; r = 0.0001;
  for (j = 0; j < n; j++){
    rand_circle(&r, &x, &y);
    fprintf(fd, "%Lf %Lf \n", x, y);
  }
  fclose(fd);
  FILE *fg = fopen("plot.gnu", "w");
  fprintf(fg, "reset \n");
  fprintf(fg, "set terminal postscript enhanced 'Helvetica' 24 \n");
  fprintf(fg, "set output 'plot.eps' \n ");
  fprintf(fg, "plot 'plot.dat' w p pt 7 ps 0.2 lc 3 notitle \n");
  fclose(fg);
  system("gnuplot plot.gnu");
  system("evince plot.eps & \n");
}

void rdm_test() {
  unsigned long long int seed;  srand(time(NULL));  seed = rand();
  void init_genrand64();  init_genrand64(seed);
  int ns = pow(10,3), nqb = 6;
  int j, k, d;
  double coh1, coh2, coh3;
  double _Complex *rrho1, *rrho2, *rrho3;
  void rdm_ginibre(int *, double _Complex *);
  void rdm_pos(int *, double _Complex *);
  void rdm_pos_sub(int *, double _Complex *);
  double coh_l1(int *, double _Complex *);
  FILE *fd = fopen("plot.dat", "w");
  for (j = 0; j < nqb; j++) {
    d = pow(2,j+1);
    rrho1 = (double _Complex *)malloc(d*d*sizeof(double _Complex));
    rrho2 = (double _Complex *)malloc(d*d*sizeof(double _Complex));
    rrho3 = (double _Complex *)malloc(d*d*sizeof(double _Complex));
    coh1 = 0.0; coh2 = 0.0; coh3 = 0.0;
    for (k = 0; k < ns; k++) {
      rdm_ginibre(&d, rrho1); coh1 += coh_l1(&d, rrho1);
      rdm_pos(&d, rrho2); coh2 += coh_l1(&d, rrho2);
      rdm_pos_sub(&d, rrho3); coh3 += coh_l1(&d, rrho3);
    }
    coh1 /= ((double) ns); coh2 /= ((double) ns); coh3 /= ((double) ns);
    printf("%i %f %f %f \n", d, coh1, coh2, coh3);
    fprintf(fd, "%i %f %f %f \n", d, coh1, coh2, coh3);
    free(rrho1); free(rrho2); free(rrho3);
  }
  fclose(fd);
  FILE *fg = fopen("plot.gnu", "w");
  fprintf(fg, "reset \n");
  fprintf(fg, "set terminal postscript enhanced 'Helvetica' 24 \n");
  fprintf(fg, "set output 'plot.eps' \n ");
  fprintf(fg, "plot [][0:1] 'plot.dat' u 1:2 w lp,'' u 1:3 w lp,'' u 1:4 w lp \n");
  fclose(fg);
  system("gnuplot plot.gnu");
  system("evince plot.eps & \n");
}

int main() {
  void rng_init(); rng_init();
  void test_rand_circle(); test_rand_circle();
  //void rdm_test(); rdm_test();
  return 0;
}


// gcc rdmg.c mat_func.c MT19937_64.c rng.c coherence.c entropy.c lapack.c gates.c rpvg.c -llapacke -lm

