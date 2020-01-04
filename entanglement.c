
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <sys/utsname.h>

void entanglement_test(){
  double p = 0.0, dp = 0.01, Ec;
  double _Complex rhop[4][4];
  void werner2qb();
  int j, k;
  int d;  d = 4;
  double concurrence();
  void ArrayDisplayC();
  FILE *fd = fopen("ent.dat", "w");
  p -= dp;
  while(p < 0.99){
    p += dp;  werner2qb(&p, rhop);
    Ec = concurrence(&d, rhop);
    printf("%f %f \n", p, Ec);
    fprintf(fd, "%f %f \n", p, Ec);
  }
  fclose(fd);
  FILE *fg = fopen("ent.gnu", "w");
  fprintf(fg, "reset \n");
  fprintf(fg, "set terminal postscript enhanced 'Helvetica' 24 \n");
  fprintf(fg, "set output 'ent.eps' \n ");
  fprintf(fg, "plot 'ent.dat' w lp \n");
  fclose(fg);
  system("gnuplot ent.gnu");
  struct utsname osf;  uname(&osf);
  long unsigned int strlen();
  if(strlen(osf.sysname) == 5){ // Linux
    system("evince ent.eps & \n");
  } else if(strlen(osf.sysname) == 6){ // Darwin (Mac)
    system("open -a skim ent.eps & \n");
  }
}

// Returns the entanglement concurrence, for two-qubit states
double concurrence(int *d, double _Complex rho[][*d]){
  void ArrayDisplayC();
  int j, ds = 2;
  double _Complex R[*d][*d], rhot[*d][*d], rhoc[*d][*d], A[*d][*d];
  double _Complex kp[4][4] = {{0,0,0,-1},{0,0,1,0},{0,1,0,0},{-1,0,0,0}};
  void matconj();  matconj(d, d, rho, rhoc);
  void matmulC();  matmulC(d, d, d, kp, rhoc, A);  matmulC(d, d, d, A, kp, rhot);
  matmulC(d, d, d, rho, rhot, R);
  double egval[*d];  double _Complex egvalR[*d];
  char jobz = 'N';
  void lapacke_zgeev();  lapacke_zgeev(&jobz, d, R, egvalR);
  //void lapacke_zheevd();  lapacke_zheevd(&jobz, d, R, egvalR);
  for(j = 0; j < (*d); j++){
    egval[j] = creal(egvalR[j]);
  }
  double maxArray1D(); double egvalMax;  egvalMax = maxArray1D(d, egval);
  double cc = 2.0*sqrt(egvalMax) - sqrt(egval[1])
              - sqrt(egval[2]) - sqrt(egval[3]) - sqrt(egval[4]);
  double conc = 0.0;  if(cc > 0.0){conc = cc;}
  return conc;
}
