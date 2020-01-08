#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>

void rdm_ginibre(int *d, long double _Complex *rdm) {
  int j,k,l;
  long double _Complex *G; 
  G = (long double _Complex *)malloc((*d)*(*d)*sizeof(long double _Complex));
  void ginibre(int *, long double _Complex *); ginibre(d, G);
  long double N2, norm_hs(int *, long double _Complex *);
  N2 = powl(norm_hs(d,G),2.0);
  for (j = 0; j < (*d); j++) {
    for (k = j; k < (*d); k++) {
      for (l = 0; l < (*d); l++) {
        *(rdm+j*(*d)+k) += (creall(*(G+j*(*d)+l)))*(creall(*(G+k*(*d)+l)));
        *(rdm+j*(*d)+k) += (cimagl(*(G+j*(*d)+l)))*(cimagl(*(G+k*(*d)+l)));
        *(rdm+j*(*d)+k) -= I*(creall(*(G+j*(*d)+l)))*(cimagl(*(G+k*(*d)+l)));
        *(rdm+j*(*d)+k) += I*(cimagl(*(G+j*(*d)+l)))*(creall(*(G+k*(*d)+l)));
      }
      *(rdm+j*(*d)+k) /= N2;
      if (j != k) {
        *(rdm+k*(*d)+j) = creall(*(rdm+j*(*d)+k)) - I*cimagl(*(rdm+j*(*d)+k));
      }
    }
  }
}

void ginibre(int *d, long double _Complex *G){
  long double grn1, grn2;
  void rng_gauss(long double *, long double *);
  int j, k;
  for (j = 0; j < (*d); j++) {
    for (k = 0; k < (*d); k++) {
      rng_gauss(&grn1, &grn2);
      *(G+j*(*d)+k) = grn1 + I*grn2;
    }
  }
}

void rdm_pos_sub(int *d, long double _Complex *rrho) {
  long double *rpv; 
  rpv = (long double *)malloc((*d)*sizeof(long double));
  void rpv_zhsl(int *, long double *); rpv_zhsl(d, rpv);
  void rand_circle(long double *, long double *, long double *);
  long double r, x, y;
  int j, k;
  for (j = 0; j < (*d); j++) {
    *(rrho+j*(*d)+j) = *(rpv+j);
  }
  free(rpv);
  for (j = 0; j < ((*d)-1); j++) {
    for (k = (j+1); k < (*d); k++) {
      r = sqrtl(cabsl(*(rrho+j*(*d)+j))*cabsl(*(rrho+k*(*d)+k)));
      rand_circle(&r, &x, &y);
      *(rrho+j*(*d)+k) = x + I*y; *(rrho+k*(*d)+j) = x - I*y;
    }
  }
}

void rdm_pos(int *d, long double _Complex *rrho) {
  long double *rpv; 
  rpv = (long double *)malloc((*d)*sizeof(long double));
  void rpv_zhsl(int *, long double *); rpv_zhsl(d, rpv);
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
  r = sqrtl(r);
  for (j = 0; j < ((*d)-1); j++) {
    for (k = j+1; k < (*d); k++) {
      rand_circle(&r, &x, &y);
      *(rrho+j*(*d)+k) = x + I*y;
      *(rrho+k*(*d)+j) = x - I*y;
      r -= (powl(x,2.0)+powl(y,2.0));
    }
  }
}

void rand_circle(long double *r, long double *x, long double *y) {
  // returns a random point (x,y) uniformly distributed
  // within a circle of radius r
  long double th, ph, rh;
  double genrand64_real1();
  th = 2.0*M_PI*((long double) genrand64_real1());
  rh = (*r)*sqrtl(((long double)genrand64_real1()));
  (*x) = rh*cosl(th); 
  (*y) = rh*sinl(th);
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
  void plot(); plot();
}

void rdm_test() {
  void rng_init(); rng_init();
  int ns = pow(10,3), nqb = 6;
  int j, k, d;
  long double coh1, coh2, coh3;
  long double _Complex *rrho1, *rrho2, *rrho3;
  void rdm_ginibre(int *, long double _Complex *);
  void rdm_pos(int *, long double _Complex *);
  void rdm_pos_sub(int *, long double _Complex *);
  long double coh_l1(int *, long double _Complex *);
  FILE *fd = fopen("plot.dat", "w");
  for (j = 0; j < nqb; j++) {
    d = pow(2,j+1);
    rrho1 = (long double _Complex *)malloc(d*d*sizeof(long double _Complex));
    rrho2 = (long double _Complex *)malloc(d*d*sizeof(long double _Complex));
    rrho3 = (long double _Complex *)malloc(d*d*sizeof(long double _Complex));
    coh1 = 0.0; coh2 = 0.0; coh3 = 0.0;
    for (k = 0; k < ns; k++) {
      rdm_ginibre(&d, rrho1); coh1 += coh_l1(&d, rrho1);
      rdm_pos(&d, rrho2); coh2 += coh_l1(&d, rrho2);
      rdm_pos_sub(&d, rrho3); coh3 += coh_l1(&d, rrho3);
    }
    coh1 /= ((long double) ns); coh2 /= ((long double) ns); coh3 /= ((long double) ns);
    fprintf(fd,"%i %f %f %f \n", d, ((double) coh1), ((double) coh2), ((double) coh3));
    free(rrho1); free(rrho2); free(rrho3);
  }
  fclose(fd);
  void plot(); plot();
}

void rdm_test_ineq() {
  int d, j, k, l;
  long double _Complex *rrho;
  void rdm_pos(int *, long double _Complex *);
  void rdm_pos_sub(int *, long double _Complex *);
  long double coh_hs(int *, long double _Complex *);
  long double coh, ub, neg;
  for (j = 0; j < 10; j++) {
    d = pow(2,j+1);
    rrho = (long double _Complex *)malloc(d*d*sizeof(long double _Complex));
    neg = 1.0;
    while (neg > 0.0) {
      rdm_pos(&d, rrho);
      //rdm_pos_sub(&d, rrho);
      coh = coh_hs(&d, rrho);
      ub = 0.0;
      for (k = 0; k < (d-1); k++) {
       for (l = (k+1); l < d; l++) {
          ub += (*(rrho+k*d+k))*(*(rrho+l*d+l));
        }
      }
      ub *= 2;
      neg = ub-coh;
    }
    printf("%i %f \n", d, ((double) neg));
    free(rrho);
  }
}

int main() {
  void rng_init(); rng_init();
  //void test_rand_circle(); test_rand_circle();
  //void rdm_test(); rdm_test();
  void rdm_test_ineq(); rdm_test_ineq();
  return 0;
}


// gcc rdmg.c mat_func.c MT19937_64.c rng.c coherence.c entropy.c lapack.c gates.c rpvg.c plots.c -llapacke -lm

