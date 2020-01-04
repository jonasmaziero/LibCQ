//------------------------------------------------------------------------------
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
//#include <math.h>
//------------------------------------------------------------------------------
void pTraceTest(){
  int da = 2, db = 2, dc = 2, d = da*db*dc, dac = da*dc;
  double _Complex rho[d][d];  // Arrays must be initialized in C
  int j, k;
  double _Complex rho_b[db][db], rho_a[da][da], rho_ac[da*dc][da*dc];
  int da = 2, db = 4, d = da*db;
  double _Complex rho[d][d];  // Arrays must be initialized in C
  int j, k;
  double _Complex rho_b[db][db], rho_a[da][da];

  for (j = 0; j < d; j++){
    for (k = 0; k < d; ++k){
      rho[j][k] = 0.0;
    }
  }
  rho[1][1] = 1.0/3.0;  rho[1][2] = 1.0/3.0;  rho[1][4] = 1.0/3.0;
  rho[2][1] = 1.0/3.0;  rho[2][2] = 1.0/3.0;  rho[2][4] = 1.0/3.0;
  rho[4][1] = 1.0/3.0;  rho[4][2] = 1.0/3.0;  rho[4][4] = 1.0/3.0;

  void array2DisplayC();
  //array2DisplayC(&d, &d, rho);
  //void pTraceR();  pTraceR(&da, &db, rho, rho_b);
  //array2DisplayC(&db, &db, rho_b);
  //void pTraceL();  pTraceL(&da, &db, rho, rho_a);
  //array2DisplayC(&da, &da, rho_a);
  //void pTraceI();  pTraceI(&da, &db, &dc, rho, rho_ac);
  //array2DisplayC(&dac, &dac, rho_ac);
  int nss = 2;
  int dr = db*dc;
  int di[2] = {da, dr};
  int ssys[2] = {0, 1};
  double _Complex rhor[dr][dr];
  void pTrace();  pTrace(&nss, di, ssys, rho, rhor);
  array2DisplayC(&dr, &dr, rhor);
}
//------------------------------------------------------------------------------
void pTrace(int *nss, int *di, int *ssys,
            double _Complex *rho, double _Complex *rhor) {
// Returns the partial trace for general multipartite systems
// nss  == No. of subsystems
// di(0:nss-1)  == dimensions of the sub-systems
// ssys(0:nss-1)  == (components = 0 or 1) subsystems to be traced out
//   If ssys(j) = 0 the j-th subsystem is traced out
//   If ssys(j) = 1 it is NOT traced out
// rho == Total density matrix
// rhor == Reduced matrix
int j, k, l, da, db, dc, dn;
// d == Total dimension (is the product of the subsystems dimensions)
// dr == Dimension of the reduced matrix (is the product of the dimensions
// of the subsystems we shall not trace out)
int d = 1;  // computing the total and reduced state dimensions
int dr = 1;
for (j = 0; j < (*nss); j++) {
  d *= (*(di+j));
  if ((*(ssys+j)) == 1) {
    dr *= (*(di+j));
  }
}
if (dr == d) {  // warning
  printf("problem: dr = d \n");
  return;
}
double _Complex *rhored;
size_t szc = sizeof(double _Complex);
void pTraceL();
void pTraceR();
void pTraceI();
void matA2Bc(); //(int *xd, int *yd, double _Complex *A, double _Complex *B);
// For bipartite systems
if ((*nss) == 2) {
  if ((*ssys) == 0) {
    pTraceL(di, di+1, rho, rhor);
  } else if ((*(ssys+1)) == 0) {
    pTraceR(di, di+1, rho, rhor)  ;
  }
  rhored = (double _Complex *)malloc(dr*dr*szc);
  matA2Bc(&dr, &dr, rhor, rhored);
// For multipartite systems
} else if (nss >= 3) {
  // Left partial traces
  l = -1;  // l defines up to which position we shall trace out, by the lhs
  while (ssys[l+1] == 0) { l += 1; }
  if ( l == -1 ) { // We do not take the left partial trace in this case
    dn = d ;   allocate(mat1(1:dn,1:dn))
    mat1 = rho  // This matrix shall be used below if l = 0
  } else if ( l > 0 ) {  // Taking left partial trace
    if ( l == 0 ) {
      da = di(1) ;
    } else {
      da = product(di(1:l)) ;
    }
    db = d/da
    allocate(mat1(1:db,1:db)) ;
    call partial_trace_a_he(rho, da, db, mat1)
    dn = db  // After this operation, the matrix left over is mat1,
             //whose dimension is dn = db
  }
  // Right partial traces
  k = nss+1 ;   do ;   if ( ssys(k-1) == 1 ) exit ;   k = k - 1 ;   enddo
  // k defines up to which position we shall trace out
  if ( k == (nss+1) ) then // We'll not take the right pTrace in this case
    allocate(mat2(1:dn,1:dn))
    mat2 = mat1  // This matrix shall be used below if k = nss+1
    allocate(rhored(1:dn,1:dn)) ;   rhored = mat2 ;   deallocate(mat1)
  else if ( k < (nss+1) ) then  // Taking the right partial trace
    if ( k == nss ) then ;   db = di(nss) ;   else ;   db = product(di(k:nss)) ;
    endif ;   da = d/db
    allocate( mat2(1:da,1:da) ) ;   call partial_trace_b_he(mat1, da, db, mat2)
    dn = da  // After this operation, the matrix left over is mat2,
    //whose dimension is dn = da
    allocate(rhored(1:dn,1:dn)) ;   rhored = mat2 ;   deallocate( mat1 )
  endif
// Inner partial traces
  if ( (k-l) > 3 ) then  //If (k-l)<=3 there is no need to take inner pTraces
  do j = (l+2), (k-2)
    if ( ssys(j) == 0 ) then
      deallocate(rhored)
      db = di(j) ;   if ( j == (k-2) ) then ;   dc = di(k-1) ;   else ;
      dc = product(di(j+1:k-1)) ;  endif ;   da = dn/(db*dc)
      allocate(mat1(1:da*dc,1:da*dc)) ;
      call partial_trace_3_he(mat2, da, db, dc, mat1)
      allocate(rhored(1:da*dc,1:da*dc)) ;   rhored = mat1
      if (j < (k-2)) then ;   dn = da*dc ;   deallocate(mat2) ;
      allocate(mat2(1:dn,1:dn)) ;   mat2 = mat1 ;   endif
      deallocate(mat1)
    endif
  enddo
  endif
  deallocate(mat2)
}
matA2Bc(&dr, &dr, rhored, rhor);
free(rhored);
}
//------------------------------------------------------------------------------
void pTraceR(int *da, int *db, double _Complex *AB, double _Complex *A) {
// Returns the right partial trace over B, for the computational basis
// representation of the global density opertator AB
int ja, ka, jb, jau, kau, jaau;
int d = (*da)*(*db);
for (ja = 0; ja < (*da); ja++) {
  jau = ja*(*db);
  jaau = ja*(*da);
  for (ka = ja; ka < (*da); ka++) {
    kau = ka*(*db);
    *(A+ka+jaau) = 0.0;
    for (jb = 0; jb < (*db); jb++) {
      *(A+ka+jaau) += *(AB+kau+jb+(jau+jb)*d);
    }
    if (ja != ka) {
      *(A+ja+ka*(*da)) = creal(*(A+ka+jaau)) - I*cimag(*(A+ka+jaau));
    }
  }
}
}
//------------------------------------------------------------------------------
void pTraceL(int *da, int *db, double _Complex *AB, double _Complex *B) {
// Returns the left partial trace over A, for the computational basis
// representation of the global density opertator AB
int ja, jb, kb, jbau, jaau;
int d = (*da)*(*db);
for (jb = 0; jb < (*db); jb++) {
  jbau = jb*(*db);
  for (kb = jb; kb < (*db); kb++) {
    *(B+kb+jbau) = 0.0;
    for (ja = 0; ja < (*da); ja++) {
      jaau = ja*(*db);
      *(B+kb+jbau) += *(AB+jaau+kb+(jaau+jb)*d);
    }
    if (jb != kb) {
      *(B+jb+kb*(*db)) = creal(*(B+kb+jbau)) - I*cimag(*(B+kb+jbau));
    }
  }
}
}
//------------------------------------------------------------------------------
void pTraceI(int *da, int *db, int *dc,
             double _Complex *ABC, double _Complex *AC) {
// Returns the inner partial trace over B, for the computational basis
// representation of the global density opertator ABC
int j, k, l, m, o, cj, ck, ccj, cck, kdc, jdc;
int dac = (*da)*(*dc), dbc = (*db)*(*dc), d = (*db)*dac;
for (j = 0; j < (*da); j++) {
  jdc = j*(*dc);
  for (l = 0; l < (*dc); l++) {
    cj = jdc+l;
    ccj = jdc*(*db)+l;
    for (m = 0; m < (*da); m++) {
      for (o = 0; o < (*dc); o++) {
        ck = m*(*dc)+o;
        cck = m*dbc+o;
        *(AC+ck+cj*dac) = 0.0 + I*0.0;
        for (k = 0; k < (*db); k++) {
          kdc = k*(*dc);
          *(AC+ck+cj*dac) += *(ABC+cck+kdc+(ccj+kdc)*d);
  //array_display(&d, &d, rho);
  //void partial_trace_a();  partial_trace_a(&da, &db, rho, rho_b);
  //ArrayDisplayC(&db, &db, rho_b);
  void pTraceB();  pTraceB(&da, &db, rho, rho_a);
  ArrayDisplayC(&da, &da, rho_a);
}
!-------------------------------------------------------------------------------
void pTrace(int *d, int *dr, int *nss, int *di, int *ssys,
            double _Complex rho[][*d], double _Complex rhor[][*dr]) {
  // Returns the partial trace for general multipartite systems
  // nss  == No. of subsystems
  // di(0:nss-1)  == dimensions of the sub-systems
  // ssys(0:nss-1)  == (components = 0 or 1) subsystems to be traced out
  //   If ssys(j) = 0 the j-th subsystem is traced out
  //   If ssys(j) = 1 it is NOT traced out
  // d == Total dimension (is the product of the subsystems dimensions)
  // rho == Total density matrix
  // dr == Dimension of the reduced matrix (is the product of the dimensions
  //   of the subsystems we shall not trace out)
  // rhor == Reduced matrix

  //complex(8), allocatable :: mat1(:,:), mat2(:,:), rhored(:,:)
  int j, k, l, da, db, dc, dn;

  // For bipartite systems
  if (nss == 2) {
    if (ssys[0] == 0) {
      da = di[0];  db = di[1];
      pTraceA(&da, &db, rho, rhor);

      allocate(rhored(1:di(2),1:di(2)))
    } else if (ssys(2) == 0) {
      call partial_trace_b_he(rho, di(1), di(2), rhor)  ;
      allocate(rhored(1:di(1),1:di(1)))
    }
   rhored = rhor
  // For multipartite systems
  } else if (nss >= 3) {
    // Left partial traces
    l = 0 ;   do ;   if ( ssys(l+1) == 1 ) exit ;   l = l + 1 ;   enddo
    ! l defines up to which position we shall trace out
    if ( l == 0 ) then // We do not take the left partial trace in this case
      dn = d ;   allocate(mat1(1:dn,1:dn))
      mat1 = rho  // This matrix shall be used below if l = 0
    else if ( l > 0 ) then  // Taking left partial trace
      if ( l == 1 ) then ;   da = di(1) ;   else ;   da = product(di(1:l)) ;
      endif ;   db = d/da
      allocate(mat1(1:db,1:db)) ;   call partial_trace_a_he(rho, da, db, mat1)
      dn = db  // After this operation, the matrix left over is mat1,
      //whose dimension is dn = db
    endif
    // Right partial traces
    k = nss+1 ;   do ;   if ( ssys(k-1) == 1 ) exit ;   k = k - 1 ;   enddo
    ! k defines up to which position we shall trace out
    if ( k == (nss+1) ) then // We shall not take the right partial trace in this case
      allocate(mat2(1:dn,1:dn))
      mat2 = mat1  ! This matrix shall be used below if k = nss+1
      allocate(rhored(1:dn,1:dn)) ;   rhored = mat2 ;   deallocate(mat1)
    else if ( k < (nss+1) ) then  // Taking the right partial trace
      if ( k == nss ) then ;   db = di(nss) ;   else ;   db = product(di(k:nss)) ;
      endif ;   da = d/db
      allocate( mat2(1:da,1:da) ) ;   call partial_trace_b_he(mat1, da, db, mat2)
      dn = da  // After this operation, the matrix left over is mat2,
      //whose dimension is dn = da
      allocate(rhored(1:dn,1:dn)) ;   rhored = mat2 ;   deallocate( mat1 )
    endif
  ! Inner partial traces
    if ( (k-l) > 3 ) then  // If (k-l) <= 3 there is no need to take inner partial traces
    do j = (l+2), (k-2)
      if ( ssys(j) == 0 ) then
        deallocate(rhored)
        db = di(j) ;   if ( j == (k-2) ) then ;   dc = di(k-1) ;   else ;
        dc = product(di(j+1:k-1)) ;  endif ;   da = dn/(db*dc)
        allocate(mat1(1:da*dc,1:da*dc)) ;
        call partial_trace_3_he(mat2, da, db, dc, mat1)
        allocate(rhored(1:da*dc,1:da*dc)) ;   rhored = mat1
        if (j < (k-2)) then ;   dn = da*dc ;   deallocate(mat2) ;
        allocate(mat2(1:dn,1:dn)) ;   mat2 = mat1 ;   endif
        deallocate(mat1)
      endif
    enddo
    endif
    deallocate(mat2)
  }
rhor = rhored ;   deallocate(rhored)

}
//------------------------------------------------------------------------------
void pTraceB(int *da, int *db, double _Complex rhoAB[][(*da)*(*db)],
             double _Complex rhoA[][*da]) {
  // Returns the right partial trace (over B), for the computational basis
  // representation of the global density opertator AB
  int ja, ka, jb, jaux, kaux;

  for (ja = 0; ja < (*da); ++ja) {
    jaux = ja*(*db);
    for (ka = ja; ka < (*da); ++ka) {
      kaux = ka*(*db);
      rhoA[ja][ka] = 0.0;
      for (jb = 0; jb < (*db); ++jb) {
        rhoA[ja][ka] += rhoAB[jaux+jb][kaux+jb];
      }
      if (ja != ka) rhoA[ka][ja] = creal(rhoA[ja][ka]) - I*cimag(rhoA[ja][ka]);
    }
  }

}
//------------------------------------------------------------------------------
void pTraceA(int *da, int *db, double _Complex rhoAB[][(*da)*(*db)],
             double _Complex rhoB[][*db]) {
  // Returns the left partial trace (over A), for the computational basis
  // representation of the global density opertator AB
  int ja, jb, kb;

  for (jb = 0; jb < (*db); ++jb) {
    for (kb = jb; kb < (*db); ++kb) {
      rhoB[jb][kb] = 0.0;
      for (ja = 0; ja < (*da); ++ja) {
        rhoB[jb][kb] += rhoAB[ja*(*db)+jb][ja*(*db)+kb];
      }
      if (jb != kb) rhoB[kb][jb] = creal(rhoB[jb][kb]) - I*cimag(rhoB[jb][kb]);
    }
  }

}
//------------------------------------------------------------------------------
void pTrace3(int *da, int *db, int *dc,
             double _Complex rhoABC[][(*da)*(*db)*(*dc)],
             double _Complex rhoAC[][(*da)*(*dc)]) {
  // Returns the inner partial trace (over B), for the computational basis
  // representation of the global density opertator ABC
  int j, k, l, m, o, cj, ck, ccj, cck;

  for (j = 0; j < (*da); j++) {
    for (l = 0; l < (*dc); l++) {
      cj = j*(*dc) + l;  ccj = j*(*db)*(*dc) + l;
      for (m = 0; m < (*da); m++) {
        rhoAC[cj][ck] = 0.0 + I*0.0;
        for (o = 0; o < (*dc); o++) {
          ck = m*(*dc) + o;  cck = m*(*db)*(*dc) + o;
          for (k = 0; k < (*db); k++) {
            rhoAC[cj][ck] += rhoABC[ccj + k*(*dc)][cck + k*(*dc)];
          }
        }
      }
    }
  }
}
}
//------------------------------------------------------------------------------
