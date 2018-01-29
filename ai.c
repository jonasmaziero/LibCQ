//-----------------------------------------------------------------------------------------------------------------------------------
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/utsname.h>
//-----------------------------------------------------------------------------------------------------------------------------------
void hwdigits_mnist(){ // We use the algorithm described in Nielsen's deep learning book
  int Ni = 784, Nh1 = 16, Nh2 = 16, No = 10;  // number of neurons in the input, hidden and output layers
  double eta = 3.0; // learning rate
  double dw = 0.01; // change rate for the NN parameters
  int ntr = pow(10,3), nte = pow(10,3);  // number of training and test images
  
  printf("Reading MNIST \n");
  double trI[ntr][Ni], teI[nte][Ni];  // arrays for the training and test images
  int trL[ntr], teL[nte];  // arrays for the labels of the train and test images
  void mnistRead();  mnistRead(&ntr, &nte, &Ni, trI, trL, teI, teL);
  
  printf("Initializing the NN \n");
  // total number of parameters to ajust: N = Ni*Nh1 + Nh1*Nh2 + Nh2*No + Nh1 + Nh2 + No (= 13002 in this case)
  double Wih1[Ni][Nh1], Wh1h2[Nh1][Nh2], Wh2o[Nh2][No], Bh1[Nh1], Bh2[Nh2], Bo[No];  // weights and biases
  void NNinit();  NNinit(&Ni, &Nh1, &Nh2, &No, Wih1, Wh1h2, Wh2o, Bh1, Bh2, Bo);  //void ArrayDisplayR();  ArrayDisplayR(&Nh2, &No, Who);
  
  printf("Testing the NN randomly initialized \n");
  void NNtest();
  NNtest(&Ni, &Nh1, &Nh2, &No, Wih1, Wh1h2, Wh2o, Bh1, Bh2, Bo, &nte, teI, teL);
  /*printf("Showing an especified digit on the screen, from an .eps file \n");
  int dp = 4;  // the digit one wants to see (index whithin the related interval)
  void digitShow();  digitShow(&dp, &teL[dp], &Ni, teI);*/
  
  printf("Training the NN \n");
  void NNtrain();
  NNtrain(&Ni, &Nh1, &Nh2, &No, Wih1, Wh1h2, Wh2o, Bh1, Bh2, Bo, &ntr, trI, trL, &eta, &dw);
    
  printf("Testing the NN \n");
  NNtest(&Ni, &Nh1, &Nh2, &No, Wih1, Wh1h2, Wh2o, Bh1, Bh2, Bo, &nte, teI, teL);
  
}
//-----------------------------------------------------------------------------------------------------------------------------------
// Initializes the weights and biases to Gaussian distributed random numbers 
void NNinit(int *Ni, int *Nh1, int *Nh2, int *No, 
            double Wih1[][*Nh1], double Wh1h2[][*Nh2], double Wh2o[][*No], double *Bh1, double *Bh2,  double *Bo){
  double rng_gauss();
  void rng_init();  rng_init();  // initializes the MT random number generator
  int j, k;
  for(j = 0; j < (*Ni); j++){
    for(k = 0; k < (*Nh1); k++){
      Wih1[j][k] = rng_gauss();
    }
  }
  for(j = 0; j < (*Nh1); j++){
    for(k = 0; k < (*Nh2); k++){
      Wh1h2[j][k] = rng_gauss();
    }
  }
  for(j = 0; j < (*Nh2); j++){
    for(k = 0; k < (*No); k++){
      Wh2o[j][k] = rng_gauss();
    }
  }
  for(j = 0; j < (*Nh1); j++){
    Bh1[j] = rng_gauss();
  }
  for(j = 0; j < (*Nh2); j++){
    Bh2[j] = rng_gauss();
  }
  for(j = 0; j < (*No); j++){
    Bo[j] = rng_gauss();
  }
}
//-----------------------------------------------------------------------------------------------------------------------------------
void mnistRead(int *ntr, int *nte, int *Ni, double trI[][*Ni], int *trL, double teI[][*Ni], int *teL){
  int j, k;
  // Writes and run a Python code that "reads" a given number of images from NMIST
  //FILE *fp = fopen("mnistNs", "w");
  //fprintf(fp, "%d \n", (*ntr));
  //fprintf(fp, "%d \n", (*nte));
  //fclose(fp);
  //system("python3 mnistRead.py");  //system("gedit labels_train &");
  // Opening the file written by python and storing the data into arrays
  FILE *fltr = fopen("trainL", "r");
  for(j = 0; j < (*ntr); j++){
    fscanf(fltr,"%d \n", &trL[j]);    
  }
  fclose(fltr);
  FILE *fdtr = fopen("trainD", "r");
  for(j = 0; j < (*ntr); j++){
    for(k = 0; k < (*Ni); k++){
      fscanf(fdtr,"%lf \t", &trI[j][k]);
    }
    fscanf(fdtr, "\n");
  }
  fclose(fdtr);
  FILE *flte = fopen("testL", "r");
  for(j = 0; j < (*nte); j++){
    fscanf(flte,"%d \n", &teL[j]);    
  }
  fclose(flte);
  FILE *fdte = fopen("testD", "r");
  for(j = 0; j < (*nte); j++){
    for(k = 0; k < (*Ni); k++){
      fscanf(fdte,"%lf \t", &teI[j][k]);
    }
    fscanf(fdte, "\n");
  }
  fclose(fdte);
}
//-----------------------------------------------------------------------------------------------------------------------------------
void digitShow(int *dp, int *dig, int *Ni, double txI[][*Ni]){
  // dp is the digit place (from the given list of digits) one wants to see
  // dl is the digit label
  printf("%d \n", *dig);
  FILE *fdi = fopen("digit.dat", "w");
  int k;
  for(k = 0; k < (*Ni); k++){
    fprintf(fdi, "%f \t", txI[*dp][k]);
    if (k > 0 && (k+1) % 28 == 0) fprintf(fdi, "\n");
  }
  fclose(fdi);
  // writing a gnuplot script, from C, and seeing the plot
  FILE *fg = fopen("digit.gnu", "w");
  fprintf(fg, "reset \n");
  fprintf(fg, "set terminal postscript enhanced 'Helvetica' 24 \n");
  fprintf(fg, "set output 'digit.eps' \n");
  fprintf(fg, "set pm3d map clip4in corners2col c1 \n");
  fprintf(fg, "splot 'digit.dat' matrix using 1:(1-$2):3 with pm3d notitle \n");
  fclose(fg);
  system("gnuplot digit.gnu");
  struct utsname osf;  uname(&osf);  //printf("%s \n", osf.sysname);
  long unsigned int strlen();  //printf("%lu \n", strlen(osf.sysname));
  //printf("%s \n", osf.sysname);
  if(strlen(osf.sysname) == 5){ // Linux
    system("evince digit.eps & \n");
  } else if(strlen(osf.sysname) == 6){ // Darwin (Mac)
    system("open -a skim digit.eps & \n");
  } else{ // Windows, etc
    printf("Are you using Windows? Really? Are you ok? \n");
  }
}
//-----------------------------------------------------------------------------------------------------------------------------------
void NNtrain(int *Ni, int *Nh1, int *Nh2, int *No,
             double Wih1[][*Nh1], double Wh1h2[][*Nh2], double Wh2o[][*No], double *Bh1, double *Bh2,  double *Bo,
             int *ntr, double trI[][*Ni], int *trL, double *eta, double *dw){ 
double dC, C, Cp, cost();
double Wih1N[*Ni][*Nh1], Wh1h2N[*Nh1][*Nh2], Wh2oN[*Nh2][*No], Bh1N[*Nh1], Bh2N[*Nh2], BoN[*No];
int j, k, l = 0;
while(l < 10){l++;
  for(j = 0; j < (*Ni); j++){
    for(k = 0; k < (*Nh1); k++){
      Wih1N[j][k] = Wih1[j][k];
    }
  }
  for(j = 0; j < (*Nh1); j++){
    for(k = 0; k < (*Nh2); k++){
      Wh1h2N[j][k] = Wh1h2[j][k];
    }
  }
  for(j = 0; j < (*Nh2); j++){
    for(k = 0; k < (*No); k++){
      Wh2oN[j][k] = Wh2o[j][k];
    }
  }
  for(j = 0; j < (*Nh1); j++){
    Bh1N[j] = Bh1[j];
  }
  for(j = 0; j < (*Nh2); j++){
    Bh2N[j] = Bh2[j];
  }
  for(j = 0; j < (*No); j++){
    BoN[j] = Bo[j];
  }
  C = cost(Ni, Nh1, Nh2, No, Wih1N, Wh1h2N, Wh2oN, Bh1N, Bh2N, BoN, ntr, trI, trL);
  printf("cost = %f \n", C);
  // updating the weights and biases to a new point of the manifold
  for(j = 0; j < (*Ni); j++){
    for(k = 0; k < (*Nh1); k++){
      Wih1N[j][k] += (*dw);
      Cp = cost(Ni, Nh1, Nh2, No, Wih1N, Wh1h2N, Wh2oN, Bh1N, Bh2N, BoN, ntr, trI, trL);
      dC = Cp-C;  Wih1[j][k] -= (*eta)*(dC/(*dw));
      Wih1N[j][k] -= (*dw);
    }
  }
  for(j = 0; j < (*Nh1); j++){
    for(k = 0; k < (*Nh2); k++){
      Wh1h2N[j][k] += (*dw);
      Cp = cost(Ni, Nh1, Nh2, No, Wih1N, Wh1h2N, Wh2oN, Bh1N, Bh2N, BoN, ntr, trI, trL);
      dC = Cp-C;  Wh1h2[j][k] -= (*eta)*(dC/(*dw));
      Wh1h2N[j][k] -= (*dw);
    }
  }
  for(j = 0; j < (*Nh2); j++){
    for(k = 0; k < (*No); k++){
      Wh2oN[j][k] += (*dw);
      Cp = cost(Ni, Nh1, Nh2, No, Wih1N, Wh1h2N, Wh2oN, Bh1N, Bh2N, BoN, ntr, trI, trL);
      dC = Cp-C;  Wh2o[j][k] -= (*eta)*(dC/(*dw));
      Wh2oN[j][k] -= (*dw);
    }
  }
  for(j = 0; j < (*Nh1); j++){
    Bh1N[j] += (*dw);
    Cp = cost(Ni, Nh1, Nh2, No, Wih1N, Wh1h2N, Wh2oN, Bh1N, Bh2N, BoN, ntr, trI, trL);
    dC = Cp-C;  Bh1[j] -= (*eta)*(dC/(*dw));
    Bh1N[j] -= (*dw);
  }
  for(j = 0; j < (*Nh2); j++){
    Bh2N[j] += (*dw);
    Cp = cost(Ni, Nh1, Nh2, No, Wih1N, Wh1h2N, Wh2oN, Bh1N, Bh2N, BoN, ntr, trI, trL);
    dC = Cp-C;  Bh2[j] -= (*eta)*(dC/(*dw));
    Bh2N[j] -= (*dw);
  }for(j = 0; j < (*No); j++){
    BoN[j] += (*dw);
    Cp = cost(Ni, Nh1, Nh2, No, Wih1N, Wh1h2N, Wh2oN, Bh1N, Bh2N, BoN, ntr, trI, trL);
    dC = Cp-C;  Bo[j] -= (*eta)*(dC/(*dw));
    BoN[j] -= (*dw);
  }
}
C = cost(Ni, Nh1, Nh2, No, Wih1, Wh1h2, Wh2o, Bh1, Bh2, Bo, ntr, trI, trL);
printf("cost = %f \n", C);
}
//-----------------------------------------------------------------------------------------------------------------------------------
double cost(int *Ni, int *Nh1, int *Nh2, int *No,
             double Wih1[][*Nh1], double Wh1h2[][*Nh2], double Wh2o[][*No],
             double *Bh1, double *Bh2,  double *Bo,
            int *ntr, double trI[][*Ni], int *trL){
  double C = 0.0, in[*Ni], out[*No];
  int j, k;
  void NNoutput();
  for(j = 0; j < (*ntr); j++){
    for(k = 0; k < (*Ni); k++){in[k] = trI[j][k];}
    NNoutput(Ni, Nh1, Nh2, No, Wih1, Wh1h2, Wh2o, Bh1, Bh2, Bo, in, out);
    // computing the squared distance between the input and output of the NN
    out[trL[j]] -= 1.0;  for(k = 0; k < (*No); k++){C += pow(out[k],2.0);}
  }
  C /= (2.0*(*ntr));
  return C;
}
//-----------------------------------------------------------------------------------------------------------------------------------
void NNtest(int *Ni, int *Nh1, int *Nh2, int *No, 
            double Wih1[][*Nh1], double Wh1h2[][*Nh2], double Wh2o[][*No], double *Bh1, double *Bh2, double *Bo, 
            int *nte, double teI[][*Ni], int *teL){
  double in[*Ni], out[*No], max;
  int j, k, l, im;
  int nci = 0;  // number of correct identifications
  void NNoutput();
  for(j = 0; j < (*nte); j++){
    for(k = 0; k < (*Ni); k++){in[k] = teI[j][k];}
    NNoutput(Ni, Nh1, Nh2, No, Wih1, Wh1h2, Wh2o, Bh1, Bh2, Bo, in, out);
    max = 0.0;
    for(l = 0; l < (*No); l++){
      if(out[l] > max){max = out[l];  im = l;}
    }
    if(im == teL[j]){nci += 1;}
  }
  printf("nte = %d,  nci = %d,  pci = %f\n", *nte, nci, (double)nci/(*nte));
}
//-----------------------------------------------------------------------------------------------------------------------------------
void NNoutput(int *Ni, int *Nh1, int *Nh2, int *No,
              double Wih1[][*Nh1], double Wh1h2[][*Nh2], double Wh2o[][*No], double *Bh1, double *Bh2, double *Bo,
              double *in, double *out){ 
  // given the input (in), and the weigths and biases, it returns the output of the neural network (out)
  double oh1[*Nh1], oh2[*Nh2], wh1[*Ni], wh2[*Nh1], wo[*Nh2], sigmoid(), b;
  int j, k;
  for(j = 0; j < (*Nh1); j++){
    for(k = 0; k < (*Ni); k++){
      wh1[k] = Wih1[k][j];
    }
    b = Bh1[j];  oh1[j] = sigmoid(Ni, wh1, in, &b);  //printf("%f %f %f %f \n", wh[j], in[j], b, oh[j]);
  }
  for(j = 0; j < (*Nh2); j++){
    for(k = 0; k < (*Nh1); k++){
      wh2[k] = Wh1h2[k][j];
    }
    b = Bh2[j];  oh2[j] = sigmoid(Nh1, wh2, oh1, &b);  //printf("%f %f %f %f \n", wo[j], oh[j], b, out[j]);
  }
  for(j = 0; j < (*No); j++){
    for(k = 0; k < (*Nh2); k++){
      wo[k] = Wh2o[k][j];
    }
    b = Bo[j];  out[j] = sigmoid(Nh2, wo, oh2, &b);  //printf("%f %f %f %f \n", wo[j], oh[j], b, out[j]);
  }
}
//-----------------------------------------------------------------------------------------------------------------------------------
double sigmoid(int *Np, double *w, double *x, double *b){ // is the output of a sigmoid neuron
  // Np is the number of neurons of the previous layer (input or hidden)
  double sig, zz;
  int j;
  double ip = 0.0;
  //zz = innerR(Np, w, x) + (*b);  // did not work ??? see why ...
  for(j = 0; j < (*Np); j++){
    ip += (w[j])*(x[j]);
  }
  zz = ip + (*b);
  sig = 1.0/(1.0 + exp(-zz)); //printf("%f %f \n", zz, sig);
  return sig;
}
//-----------------------------------------------------------------------------------------------------------------------------------