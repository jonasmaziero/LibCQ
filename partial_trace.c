//-----------------------------------------------------------------------------------------------------------------------------------
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
//-----------------------------------------------------------------------------------------------------------------------------------
// Returns the left partial trace (over a), for a bi-partite HERMITIAN matrix rho
int partial_trace_b(int* da, int* db, double _Complex rho[][(*da)*(*db)], double _Complex rho_a[][*da])
{
// da, db ! Dimensions of the subsystems (the dimension of the whole system is d = da*db)
// rho(da*db,da*db)  ! Bipartite matrix (computational basis representation of the ragarded operator)
// rho_a(da,da)  !  Reduced matrix
int ja, ka, jb, jaux, kaux;  // Auxiliary variables for counters  
for(ja = 0; ja < (*da); ++ja)
{
  for(ka = ja; ka < (*da); ++ka)
  {
    jaux = ja*(*db);
    kaux = ka*(*db);
    rho_a[ja][ka] = 0.0;
    for(jb = 0; jb < (*db); ++jb)
    {
      rho_a[ja][ka] += rho[jaux+jb][kaux+jb];
    }
    if (ja != ka) rho_a[ka][ja] = creal(rho_a[ja][ka]) -I*cimag(rho_a[ja][ka]);
  }
}
return 0;
}
//-----------------------------------------------------------------------------------------------------------------------------------
// Returns the left partial trace (over a), for a bi-partite HERMITIAN matrix rho
int partial_trace_a(int* da, int* db, double _Complex rho[][(*da)*(*db)], double _Complex rho_b[][*db])
{
// da, db ! Dimensions of the subsystems (the dimension of the whole system is d = da*db)
// rho(da*db,da*db)  ! Bipartite matrix (computational basis representation of the ragarded operator)
// rho_b(db,db)  !  Reduced matrix
int ja, jb, kb;  // Auxiliary variables for counters  
for(jb = 0; jb < (*db); ++jb)
{
  for(kb = jb; kb < (*db); ++kb)
  {
    rho_b[jb][kb] = 0.0;
    for(ja = 0; ja < (*da); ++ja)
    {
      rho_b[jb][kb] += rho[ja*(*db)+jb][ja*(*db)+kb];
    }
    if (jb != kb) rho_b[kb][jb] = creal(rho_b[jb][kb]) -I*cimag(rho_b[jb][kb]);
  }
}
return 0;
}
//-----------------------------------------------------------------------------------------------------------------------------------