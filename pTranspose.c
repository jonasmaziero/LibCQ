//------------------------------------------------------------------------------------------------------------------------------------
//#include <stdlib.h>
//------------------------------------------------------------------------------------------------------------------------------------
// Returns its PARTIAL TRANSPOSE with relation to system A
void pTransposeA(int* da, int* db, double _Complex rho[][(*da)*(*db)], double _Complex rhoTa[][(*da)*(*db)]){  
  int ja, ka, jb, kb;
  for (ja = 0; ja < (*da); ja++){
    for (ka = 0; ka < (*da); ka++){
      for (jb = 0; jb < (*db); jb++){
        for (kb = 0; kb < (*db); kb++){
          rhoTa[ka*(*db)+jb][ja*(*db)+kb] = rho[ja*(*db)+jb][ka*(*db)+kb];
        }
      }
    }
  }
}
//------------------------------------------------------------------------------------------------------------------------------------
// Returns its PARTIAL TRANSPOSE with relation to system A
void pTransposeB(int* da, int* db, double _Complex rho[][(*da)*(*db)], double _Complex rhoTb[][(*da)*(*db)]){
  int ja, ka, jb, kb;
  for (ja = 0; ja < (*da); ja++){
    for (ka = 0; ka < (*da); ka++){
      for (jb = 0; jb < (*db); jb++){
        for (kb = 0; kb < (*db); kb++){
          rhoTb[ja*(*db)+kb][ka*(*db)+jb] = rho[ja*(*db)+jb][ka*(*db)+kb];
        }
      }
    }
  }
}
//------------------------------------------------------------------------------------------------------------------------------------