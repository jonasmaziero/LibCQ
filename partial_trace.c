//-----------------------------------------------------------------------------------------------------------------------------------
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
//-----------------------------------------------------------------------------------------------------------------------------------
// Returns the left partial trace (over a), for a bi-partite HERMITIAN matrix rho
int partial_trace_a(int* da, int* db, double _Complex rho[][(*da)*(*db)], double _Complex rho_b[][*db]){
// da, db ! Dimensions of the subsystems (the dimension of the whole system is d = da*db)
// rho(da*db,da*db)  ! Bipartite matrix (computational basis representation of the ragarded operator)
// rho_b(db,db)  !  Reduced matrix
int j, k, l;  // Auxiliary variables for counters  
for(j = 0; j < (*db); ++j)
{
  for(k = j; k < (*db); ++k)
  {
    rho_b[j][k] = 0.0;
    for(l = 0; l < (*da); ++l)
    {
      rho_b[j][k] += rho[(*db)*l+j][(*db)*l+k];
    }
    if (j != k) rho_b[k][j] = creal(rho_b[j][k]) -I*cimag(rho_b[j][k]);
  }
}
  return 0;}
//-----------------------------------------------------------------------------------------------------------------------------------