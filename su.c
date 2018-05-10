//------------------------------------------------------------------------------
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
//------------------------------------------------------------------------------
int main()
{
  double _Complex *sj, *sjsj;
  size_t szdc = sizeof(double _Complex);
  sj = (double _Complex *)malloc(4*szdc);
  sjsj = (double _Complex *)malloc(4*4*szdc);
  void pauli();  pauli(2, sj);
  int d = 2, D = 4;
  void array2DisplayC();  array2DisplayC(&d, &d, sj);
  void kronecker() ;
  kronecker(&d, &d, &d, &d, sj, sj, sjsj);
  array2DisplayC(&D, &D, sjsj);
}
//------------------------------------------------------------------------------
void pauli(int j, double _Complex sj[][2])
{
  if(j == 0){
    sj[0][0] = 1;  sj[0][1] = 0;  sj[1][0] = 0;  sj[1][1] = 1;
  } else if(j == 1){
    sj[0][0] = 0;  sj[0][1] = 1;  sj[1][0] = 1;  sj[1][1] = 0;
  } else if(j == 2){
    sj[0][0] = 0;  sj[0][1] = -1j;  sj[1][0] = 1j;  sj[1][1] = 0;
  } else if(j == 3){
    sj[0][0] = 1;  sj[0][1] = 0;  sj[1][0] = 0;  sj[1][1] = -1;
  }
}
//------------------------------------------------------------------------------
