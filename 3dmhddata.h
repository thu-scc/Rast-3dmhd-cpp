#ifndef _3DMHD_DATA_H_
#define _3DMHD_DATA_H_

#include "3dmhdparam.h"

#include <mpi.h>

typedef double data[nz][ny][nx];

extern int ipar[32];
extern double par[64];
/*BIG*/
extern data RU, RV, RW, RO, TT;
extern data UU, VV, WW;
extern data FU, FV, FW, FR, FT;
extern data ZRU, ZRV, ZRW, ZRO, ZTT;
extern data WW1, WW2, WW3;
extern data BX, BY, BZ;
extern data ZBX, ZBY, ZBZ;
/*AJACOBI*/
extern double EXX[nx], dxxdx[nx], d2xxdx2[nx], DDX[nx];
extern double WYY[ny], dyydy[ny], d2yydy2[ny], DDY[ny];
extern double ZEE[nz], dzzdz[nz], d2zzdz2[nz], DDZ[nz];
/*ITER*/
extern int nTotal, nStep0, nIt;
/*COMMUN*/
extern int myPE, myPEy, myPEz;
extern MPI_Datatype MPISIZE;
/*CPAR*/
extern double CV, OCV, ORE, RE, REPR, THETA, GRAV, AMPT, GAMMA;
/*CROT*/
extern double OMX, OMZ;
/*CMAG*/
extern double ORM, RM, OBETA, AMPB, BFH, BZP;
/*CPEN*/
extern double PZP, SIGMA, RKAPST, TB, RKAPM;
extern double RKAPA[nz], DKAPA[nz];
/*CPER*/
extern double TP, XP, YP, ZP, TC, QFH, HH;
/*BOUNDS*/
extern double xMax, yMax, zMax;
/*GRID*/
extern double DD, hx, h2x, hy, h2y, hz, h2z, C13, C23, C43;
/*CTIM*/
extern double DT, TIMT, TIMC, TIMI;
/*TRACE*/
extern double UMACH;
/*RUNGKU*/
extern double GAM1, GAM2, GAM3, ZETA1, ZETA2;
/*SPLINEX*/
extern double HHX[nPx + 3], SIGX[nPx], AAX[nPx + 3], BBX[nPx + 3], XPINV[nPx], xHH;
extern int iSEGx, kLOx[nPx + 3], kHIx[nPx + 3];
/*SPLINEY*/
extern double HHY[nRy + 3], SIGY[nRy], AAY[nRy + 3], BBY[nRy + 3], YPINV[nRy], yHH;
extern int iSEGy, kLOy[nRy + 3], kHIy[nRy + 3];
/*RELAX*/
extern double TSTART, TOFF, RLAX;
/*SPECIALBOUND*/
extern double TU, DZTB, DZTU;
/*FFCOM*/
extern double A, C;

#endif
