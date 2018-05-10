#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <math.h>

int main(){
  int j, k;
  int d = 4;
  double *rhoR, *rhoI;
  double _Complex *rho;
  double concurrence();
  rhoR = (double *)malloc(16*sizeof(double));
  rhoI = (double *)malloc(16*sizeof(double));
  rho = (double _Complex *)malloc(16*sizeof(double _Complex));
  FILE *flR = fopen("rhoR", "r");
  FILE *flI = fopen("rhoI", "r");
  for(j = 0; j < 4; j++){
    for(k = 0; k < 4; k++){
      fscanf(flR,"%f \t", *(rhoR+j*4+k));
      fscanf(flI,"%f \t", *(rhoI+j*4+k));
      *(rho+j*4+k) = *(rhoR+j*4+k) + (*(rhoI+j*4+k))*I;
    }
    fscanf(flR,"\n");
    fscanf(flI,"\n");
  }
  fclose(flR);
  fclose(flI);
  printf("%lf \n", concurrence(&d, rho));
  /*for(j = 0; j < 4; j++){
    for(k = 0; k < 4; k++){
      printf("%lf \t", *(rhoR+j*4+k));
    }
    printf("\n");
  }*/
  return 0;
}
