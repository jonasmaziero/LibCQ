//-----------------------------------------------------------------------------------------------------------------------------------
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
//-----------------------------------------------------------------------------------------------------------------------------------
// Returns the right partial trace (over b), for a bi-partite HERMITIAN matrix rho=rho_ab
void partial_trace_b(int* da, int* db, double _Complex rho[][(*da)*(*db)], double _Complex rho_a[][*da])
{
  int ja, ka, jb, jaux, kaux;
  for (ja = 0; ja < (*da); ++ja)
  {
    jaux = ja*(*db);
    for (ka = ja; ka < (*da); ++ka)
    {
      kaux = ka*(*db);
      rho_a[ja][ka] = 0.0;
      for (jb = 0; jb < (*db); ++jb)
      {
        rho_a[ja][ka] += rho[jaux+jb][kaux+jb];
      }
      if (ja != ka) rho_a[ka][ja] = creal(rho_a[ja][ka]) -I*cimag(rho_a[ja][ka]);
    }
  }
}
//-----------------------------------------------------------------------------------------------------------------------------------
// Returns the left partial trace (over a), for a bi-partite HERMITIAN matrix rho=rho_ab
void partial_trace_a(int* da, int* db, double _Complex rho[][(*da)*(*db)], double _Complex rho_b[][*db])
{
  int ja, jb, kb;
  for (jb = 0; jb < (*db); ++jb)
  {
    for (kb = jb; kb < (*db); ++kb)
    {
      rho_b[jb][kb] = 0.0;
      for (ja = 0; ja < (*da); ++ja)
      {
        rho_b[jb][kb] += rho[ja*(*db)+jb][ja*(*db)+kb];
      }
      if (jb != kb) rho_b[kb][jb] = creal(rho_b[jb][kb]) -I*cimag(rho_b[jb][kb]);
    }
  }
}
//-----------------------------------------------------------------------------------------------------------------------------------