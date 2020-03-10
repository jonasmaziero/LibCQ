
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>


void rdm_ginibre(int *d, double _Complex *rdm) {
  void zero_mat_c(int *, int *, double _Complex *); 
  zero_mat_c(d,d,rdm);
  int j, k, l;
  double _Complex *G; 
  G = (double _Complex *)malloc((*d)*(*d)*sizeof(double _Complex));
  void ginibre(int *, double _Complex *); 
  ginibre(d, G);
  double N2, norm_hs(int *, double _Complex *); 
  N2 = pow(norm_hs(d,G),2);
  for (j = 0; j < (*d); j++) {
    for (k = j; k < (*d); k++) {
      for (l = 0; l < (*d); l++) {
        *(rdm+j*(*d)+k) += conj(*(G+l*(*d)+j))*(*(G+l*(*d)+k));
      }
      *(rdm+j*(*d)+k) /= N2;
      if (j != k) {
        *(rdm+k*(*d)+j) = creal(*(rdm+j*(*d)+k)) - I*cimag(*(rdm+j*(*d)+k));
      }
    }
  }
  free(G);
}


void rdm_ginibre_classes(int *d, double _Complex *rdm, int *nulle) {
  // nulle: array of dimension dxd. nulle[j,k]=0 if
  // rdm[j,k] must be zero (k>j). nulle[j][k]=1 otherwise.
  void zero_mat_c(int *, int *, double _Complex *); zero_mat_c(d,d,rdm);
  int j, k, l;
  double _Complex *G; G = (double _Complex *)malloc((*d)*(*d)*sizeof(double _Complex));
  double _Complex *v; v = (double _Complex *)malloc((*d)*sizeof(double _Complex));
  double _Complex *w; w = (double _Complex *)malloc((*d)*sizeof(double _Complex));
  double _Complex ip_c(int *, double _Complex *, double _Complex *), ipc;
  double norm_c(int *, double _Complex *), vnorm2, norm2, norm_hs(int *, double _Complex *);
  void ginibre(int *, double _Complex *); ginibre(d, G);
  void array_display_c(int *, int *, double _Complex *);
  double vnorm2_c(int *, double _Complex *);
  for (j = 0; j < ((*d)-1); j++) {
    for (k = (j+1); k < (*d); k++) { 
      if (*(nulle+j*(*d)+k) == 0) {
        for (l = 0; l < (*d); l++) {
          *(v+l) = *(G+l*(*d)+j); // |C_j(G)>
          *(w+l) = *(G+l*(*d)+k); // |C_k(G)>
        }
        ipc = ip_c(d, v, w);
        vnorm2 = vnorm2_c(d, v);
        for (l = 0; l < (*d); l++) {
          *(G+l*(*d)+k) -= (ipc*(*(G+l*(*d)+j)))/vnorm2;
        }
      }
    }
  }
  norm2 = pow(norm_hs(d,G),2);
  for (j = 0; j < (*d); j++) {
    for (k = j; k < (*d); k++) {
      for (l = 0; l < (*d); l++) {
        *(rdm+j*(*d)+k) += conj(*(G+l*(*d)+j))*(*(G+l*(*d)+k));
      }
      *(rdm+j*(*d)+k) /= norm2;
      if (j != k) {
        *(rdm+k*(*d)+j) = creal(*(rdm+j*(*d)+k)) - cimag(*(rdm+j*(*d)+k))*I;
      }
    }
  }
  free(G); free(v); free(w);
}

void rdm_ginibre_classesl(int *d, double _Complex *rdm, int *nulle) {
  // nulle: array of dimension dxd. nulle[j,k]=0 if
  // rdm[j,k] must be zero (k>j). nulle[j][k]=1 otherwise.
  void zero_mat_c(int *, int *, double _Complex *); zero_mat_c(d,d,rdm);
  int j, k, l;
  long double _Complex *G; 
  G = (long double _Complex *)malloc((*d)*(*d)*sizeof(long double _Complex));
  long double _Complex *v;
  v = (long double _Complex *)malloc((*d)*sizeof(long double _Complex));
  long double _Complex *w; 
  w = (long double _Complex *)malloc((*d)*sizeof(long double _Complex));
  long double _Complex ip_cl(int *, long double _Complex *, long double _Complex *), ipc;
  long double norm_cl(int *, long double _Complex *), vnorm2;
  double norm_hsl(int *, long double _Complex *);
  double norm2;
  void ginibrel(int *, long double _Complex *); ginibrel(d, G);
  void array_display_c(int *, int *, double _Complex *);
  double vnorm2_cl(int *, long double _Complex *);
  for (j = 0; j < ((*d)-1); j++) {
    for (k = (j+1); k < (*d); k++) { 
      if (*(nulle+j*(*d)+k) == 0) {
        for (l = 0; l < (*d); l++) {
          *(v+l) = *(G+l*(*d)+j); // |C_j(G)>
          *(w+l) = *(G+l*(*d)+k); // |C_k(G)>
        }
        ipc = ip_cl(d, v, w);
        vnorm2 = vnorm2_cl(d, v);
        for (l = 0; l < (*d); l++) {
          *(G+l*(*d)+k) -= (ipc*(*(G+l*(*d)+j)))/vnorm2;
        }
      }
    }
  }
  norm2 = ((double) powl(norm_hsl(d,G),2));
  for (j = 0; j < (*d); j++) {
    for (k = j; k < (*d); k++) {
      for (l = 0; l < (*d); l++) {
        *(rdm+j*(*d)+k) += ((double _Complex) conjl(*(G+l*(*d)+j))*(*(G+l*(*d)+k)));
      }
      *(rdm+j*(*d)+k) /= norm2;
      if (j != k) {
        *(rdm+k*(*d)+j) = creal(*(rdm+j*(*d)+k)) - cimag(*(rdm+j*(*d)+k))*I;
      }
    }
  }
  free(G); free(v); free(w);
}


void ginibre(int *d, double _Complex *G) {
  double grn1, grn2;
  void rng_gauss(double *, double *);
  int j, k;
  for (j = 0; j < (*d); j++) {
    for (k = 0; k < (*d); k++) {
      rng_gauss(&grn1, &grn2);
      *(G+j*(*d)+k) = grn1 + grn2*I;
    }
  }
}

void ginibrel(int *d, long double _Complex *G) {
  long double grn1, grn2;
  void rng_gaussl(long double *, long double *);
  int j, k;
  for (j = 0; j < (*d); j++) {
    for (k = 0; k < (*d); k++) {
      rng_gaussl(&grn1, &grn2);
      *(G+j*(*d)+k) = grn1 + grn2*I;
    }
  }
}

/*
void rdm_pos_sub(int *d, double _Complex *rrho) { // only necessary condition
  double *rpv; 
  rpv = (double *)malloc((*d)*sizeof(double));
  void rpv_zhsl(int *, double *); rpv_zhsl(d, rpv);
  void rand_circle(double *, double *, double *);
  double r, x, y;
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


void rdm_pos(int *d, double _Complex *rrho) { // only necessary condition
  double *rpv; 
  rpv = (double *)malloc((*d)*sizeof(double));
  void rpv_zhsl(int *, double *); rpv_zhsl(d, rpv);
  void rand_circle(double *, double *, double *);
  double r, x, y;
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
  r = sqrt(r); // the initial radius
  for (j = 0; j < ((*d)-1); j++) {
    for (k = j+1; k < (*d); k++) {
      rand_circle(&r, &x, &y);
      *(rrho+j*(*d)+k) = x + I*y; 
      *(rrho+k*(*d)+j) = x - I*y;
      r = sqrtl(powl(r,2)-powl(x,2)-powl(y,2)); // update of the radius
    }
  }
}


void rdm_posl(int *d, double _Complex *rrho) { // only necessary condition
  void zero_mat_c(int *, int *, double _Complex *); zero_mat_c(d,d,rrho);
  void rand_circlel(long double *, long double *, long double *);
  long double r, x, y;
  long double *rpv; 
  rpv = (long double *)malloc((*d)*sizeof(long double));
  void rpv_zhsll(int *, long double *); rpv_zhsll(d, rpv);
  int j, k;
  for (j = 0; j < (*d); j++) {
    *(rrho+j*(*d)+j) = ((double _Complex) *(rpv+j));
    //printf("%.3Lf \t", *(rpv+j));
  }
  //printf("\n");
  r = 0.0;
  for (j = 0; j < ((*d)-1); j++) {
    for (k = (j+1); k < (*d); k++) {
      r += (*(rpv+j))*(*(rpv+k));
    }
  }
  r = sqrtl(r); // the initial radius
  //printf("%.20Lf \n", r);
  free(rpv);
  for (j = 0; j < ((*d)-1); j++) {
    for (k = j+1; k < (*d); k++) {
      rand_circlel(&r, &x, &y);
      *(rrho+j*(*d)+k) = ((double _Complex) (x + I*y)); 
      *(rrho+k*(*d)+j) = creal(*(rrho+j*(*d)+k)) -I*cimag(*(rrho+j*(*d)+k));
      r = sqrtl(powl(r,2)-powl(x,2)-powl(y,2)); // update of the radius
      //printf("r = %.20Lf, x = %.20Lf, y = %.20Lf \n", r, x, y);
    }
  }
}


void rand_circle(double *r, double *x, double *y) {
  // returns a random point (x,y) uniformly distributed
  // within a circle of radius r
  double th, ph, rh;
  double genrand64_real1();
  th = 2.0*M_PI*genrand64_real1();
  rh = (*r)*sqrt(genrand64_real1());
  (*x) = rh*cos(th); 
  (*y) = rh*sin(th);
}


void rand_circlel(long double *r, long double *x, long double *y) { // ok
  long double th, rh;
  long double genrand64_real1l();
  const long double pil = 3.141592653589793238;
  th = 2.0*pil*genrand64_real1l();
  rh = (*r)*sqrtl(genrand64_real1l());
  (*x) = rh*cosl(th); 
  (*y) = rh*sinl(th);
}


void test_rand_circle() {
  int n, j;
  //double r, x, y;
  long double r, x, y;
  void rand_circle(double *, double *, double *);
  void rand_circlel(long double *, long double *, long double *);
  FILE *fd = fopen("plot.dat", "w");
  n = pow(10,4); r = pow(10,-18);
  for (j = 0; j < n; j++){
    //rand_circle(&r, &x, &y);
    rand_circlel(&r, &x, &y);
    //fprintf(fd, "%.20f %.20f \n", x, y);
    fprintf(fd, "%.20Lf %.20Lf \n", x, y);
  }
  fclose(fd);
  void plot(); plot();
}


void rdm_test_ineq() {
  int d, j, k, l;
  double _Complex *rrho;
  void rdm_posl(int *, double _Complex *);
  void rdm_pos_sub(int *, double _Complex *);
  double coh_hs_(int *, double _Complex *);
  double coh, ub, neg, tr;
  clock_t t1, t2;
  double dt;
  //double _Complex trace_c(int *, double _Complex *);
  for (j = 0; j < 8; j++) {
    d = pow(2,j+1);
    rrho = (double _Complex *)malloc(d*d*sizeof(double _Complex));
    neg = 1.0;
    t1 = clock();
    dt = 0;
    while (neg > -0.00000001 && dt < 10) {
      rdm_posl(&d, rrho);
      //tr = ((double) trace_c(&d, rrho));
      //if (tr != 1.0) printf("d = %i, Tr(rrho) = %f \n", d, tr);
      //rdm_pos_sub(&d, rrho);
      coh = coh_hs_(&d, rrho);
      ub = 0.0;
      for (k = 0; k < (d-1); k++) {
       for (l = (k+1); l < d; l++) {
          ub += ((double) *(rrho+k*d+k))*((double) *(rrho+l*d+l));
        }
      }
      ub = 2*ub; neg = ub-coh;
      t2 = clock();
      dt = (double)((t2-t1)/CLOCKS_PER_SEC);
    }
    if (neg < -0.00000001) {
      printf("%i %.8f \n", d, neg);
    }
    else {
      printf("d = %i, Did not violate C_hs < up in %f seconds \n", d, dt);
    }
    free(rrho);
  }
}


void test_positivity() {
  void array_display_c(int *, int *, double _Complex *);
  void lapacke_zheevd(char *, int *, double _Complex *, double *);
  double *W;
  char jobz = 'N';
  int j, d, l;
  double veccsum(int *, double *), sevals;
  double _Complex trace_c(int *, double _Complex *), tr;
  double _Complex *rrho;
  //void rdm_posl(int *, double _Complex *);
  void rdm_ginibre_classes(int *, double _Complex *, int *);
  void rdm_ginibre(int *, double _Complex *);
  int *nulle;
  void zero_mat_i(int *, int *, int *); 
  for (j = 0; j < 8; j++) {
    d = pow(2,j+1);
    rrho = (double _Complex *)malloc(d*d*sizeof(double _Complex)); 
    //rdm_posl(&d, rrho);
    nulle = (int *)malloc(d*d*sizeof(int)); 
    zero_mat_i(&d,&d,nulle);
    rdm_ginibre_classes(&d, rrho, nulle);
    //rdm_ginibre(&d, rrho);
    W = (double *)malloc(d*sizeof(double)); lapacke_zheevd(&jobz, &d, rrho, W);
    //printf("rrho \n"); array_display_c(&d,&d,rrho);
    printf("Negative eigenvalues for d=%i \n", d);
    for (l = 0; l < d; l++) {
      if (*(W+l) < 0.0) printf("%f \t", *(W+l));
      //printf("%f \t",*(W+l));
    }
    printf("\n");
    sevals = veccsum(&d, W); tr =  trace_c(&d, rrho);
    printf("d=%i, sum evals = %f, tr = %f +i %f \n", d, sevals, creal(tr), cimag(tr));
    free(W); free(rrho); free(nulle);
  }
}
*/

void rdm_test() {
  void rng_init(); rng_init();
  int ns = pow(10,3), nqb = 6;
  int j, k, l, m, n, d;
  double coh1, coh2;
  double _Complex *rrho1, *rrho2;
  void rdm_ginibre(int *, double _Complex *);
  void rdm_ginibre_classes(int *, double _Complex *, int *);
  double coh_l1(int *, double _Complex *);
  int *nulle;
  void one_mat_i(int *, int *, int *); 
  void array_display_c(int *, int *, double _Complex *);
  int rn_int_ab(int *, int *), ct, a, b;
  void array_display_i(int *, int *, int *);
  void rdm_ginibre_classesl(int *, double _Complex *, int *);
  FILE *fd = fopen("plot.dat", "w");
  for (j = 0; j < nqb; j++) {
    d = pow(2,j+1);
    rrho1 = (double _Complex *)malloc(d*d*sizeof(double _Complex));
    rrho2 = (double _Complex *)malloc(d*d*sizeof(double _Complex));
    nulle = (int *)malloc(d*d*sizeof(int)); one_mat_i(&d,&d,nulle);
    ct = 0; a = 0; b = d-1;
    while (ct < d*(d-1)/4) { // chooses randomly the null coherences
      m = rn_int_ab(&a,&b); n = rn_int_ab(&a,&b);
      if (m != n) {
        *(nulle+m*d+n) = 0; *(nulle+n*d+m) = 0; ct += 1;
      }
    }
    printf("nulle \n"); if (d == 4) array_display_i(&d,&d,nulle);
    coh1 = 0; coh2 = 0;
    for (k = 0; k < ns; k++) {
      rdm_ginibre(&d, rrho1); coh1 += coh_l1(&d, rrho1);
      //rdm_ginibre_classes(&d, rrho2, nulle); coh2 += coh_l1(&d, rrho2);
      rdm_ginibre_classesl(&d, rrho2, nulle); coh2 += coh_l1(&d, rrho2);
    }
    printf("rrho \n"); if (d == 4) array_display_c(&d,&d,rrho2);
    coh1 /= ((double) ns); 
    coh2 /= ((double) ns); 
    fprintf(fd,"%i %f %f \n", d, coh1, coh2);
    printf("d = %i, C1 = %f, C2 = %f \n", d, coh1, coh2);
    free(nulle); free(rrho1); free(rrho2);
  }
  fclose(fd);
  void plot(); plot();
}


int main() {
  void rng_init(); rng_init();
  //void test_rand_circle(); test_rand_circle(); // ok
  //void rdm_test_ineq(); rdm_test_ineq(); // ok
  //void test_positivity(); test_positivity(); // not ok
  void rdm_test(); rdm_test();
  return 0;
}


// gcc rdmg.c mat_func.c MT19937_64.c rng.c coherence.c entropy.c lapack.c gates.c rpvg.c plots.c -llapacke -lm

