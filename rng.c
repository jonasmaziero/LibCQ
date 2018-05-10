//------------------------------------------------------------------------------
#include <stdlib.h>
#include <math.h>
#include <time.h>
//------------------------------------------------------------------------------
void rngTest(){
  unsigned long long int seed;  srand(time(NULL));  seed = rand();
  void init_genrand64();  init_genrand64(seed);
  /*int j;
  FILE *fp;  char output[] = "output.txt";  fp = fopen(output, "w+");
  double genrand64_real1(), a, b;
  int ns = pow(10,4);  //printf("%d \n", ns);
  for (j = 0; j < ns; j++){
    a = genrand64_real1();  b = genrand64_real1();
    fprintf(fp, "%f \t %f \n", a, b);  //printf("%f \t %f \n", a, b);
  }
  for (j = 0; j < 10; j++){
    fscanf(fp, "%f \t %f \n", genrand64_real1(), genrand64_real1());
    printf("%f \t %f \n", genrand64_real1(), genrand64_real1());
  }
  fclose(fp);*/

  // testing the Gaussians
  printf("Computing the probability distribution");
  // compile with: gcc tests.c rng.c MT19937_64.c -lm
  int j, l;
  int ns = pow(10,6), ni = 100;
  double xmin = -5.0, xmax = 5.0;
  double delta = (xmax-xmin)/ni;
  int ct[ni];
  for(j = 0; j < ni; j++){
    ct[j] = 0;
  }
  double rng_gauss(), grn;
  for(j = 0; j < ns; j++){
    grn = rng_gauss();
    for(l = 0; l < (ni-1); l++){
      if((grn >= (xmin + l*delta)) && (grn < (xmin + (l+1)*delta))) {
        ct[l] = ct[l] + 1;
      }
    }
  }
  FILE *fd = fopen("rn.dat", "wr");
  for(l = 0; l < (ni-1); l++){
    fprintf(fd, "%f %f \n", (xmin +(l+0.5)*delta), (double) ct[l]/ns);
  }
  fclose(fd);
  // writing a gnuplot script, from C
  FILE *fgnu = fopen("rn.gnu", "wr");
  fprintf(fgnu, "reset \n");
  fprintf(fgnu, "set terminal postscript enhanced 'Helvetica' 24 \n");
  fprintf(fgnu, "set output 'rn.eps' \n");
<<<<<<< HEAD
  fprintf(fgnu, "plot 'rn.dat',
          (1.0/sqrt(2.0*pi*(1.0**2)))*exp(-((x-0)**2)/(2.0*(1.0**2))) \n");
=======
  fprintf(fgnu, "plot 'rn.dat', (1.0/sqrt(2.0*pi*(1.0**2)))*exp(-((x-0)**2)/(2.0*(1.0**2))) \n");
>>>>>>> dadec15c06da70b95046f6c7311b09a8024f1e31
  fclose(fgnu);
  // running a gnuplot script, from C
  system("gnuplot rn.gnu");
  // seeing the result
  system("evince rn.eps&");
}
<<<<<<< HEAD
//------------------------------------------------------------------------------
void rng_init(){ // initializes the MT random number generator
  unsigned long long int seed;  void srand();  srand(time(NULL));
  int rand();  seed = rand();
  void init_genrand64();  init_genrand64(seed);
}
//------------------------------------------------------------------------------
double rng_gauss()
{
  // returns random numbers with Gaussian probability distribution
=======
//-----------------------------------------------------------------------------------------------------------------------------------
void rng_init(){ // initializes the MT random number generator
  unsigned long long int seed;  void srand();  srand(time(NULL));  int rand();  seed = rand();
  void init_genrand64();  init_genrand64(seed);
}
//-----------------------------------------------------------------------------------------------------------------------------------
double rng_gauss(){ // returns random numbers with Gaussian probability distribution
>>>>>>> dadec15c06da70b95046f6c7311b09a8024f1e31
  double rn1, rn2, grn;
  double pi = 4.0*atan(1.0);
  double logterm, angle;
  double genrand64_real1();
  rn1 = genrand64_real1();  rn2 = genrand64_real1();
<<<<<<< HEAD
  if (rn2 < 1.e-15) rn2 = 1.e-15;
  logterm = sqrt(-2.0*log(rn2));  angle = 2.0*pi*rn1;
  grn = logterm*cos(angle);
  return grn;
}
//------------------------------------------------------------------------------
=======
  if (rn2 < 1.e-15) rn2 = 1.e-15;  logterm = sqrt(-2.0*log(rn2));  angle = 2.0*pi*rn1;
  grn = logterm*cos(angle);
  return grn;
}
//-----------------------------------------------------------------------------------------------------------------------------------
>>>>>>> dadec15c06da70b95046f6c7311b09a8024f1e31
