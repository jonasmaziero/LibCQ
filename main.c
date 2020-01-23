#include <stdio.h>

int main()
{
  printf("%li %li %li \n", sizeof(float), sizeof(double), sizeof(long double));
  printf("%li %li %li \n", sizeof(_Complex), sizeof(double _Complex), sizeof(long double _Complex));
  return 0;
}
