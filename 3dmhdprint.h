#ifndef _3DMHD_PRINT_H_
#define _3DMHD_PRINT_H_

#include "3dmhddata.h"
#include "3dmhdparam.h"

FILE* fp;

#define for_i(n)  for (int i = 0; i < n; i++)
#define printi(a) fprintf(fp, "%d ", a)
#define printd(a) fprintf(fp, "%11.4E", a)
#define endl()    fputs("", fp)

// extern int ipar[32];
// extern double par[64];
void print_par() {
    for_i(32) printi(ipar[i]);
    for_i(64) printd(par[i]);
    endl();
}

/*BIG*/
// extern data RU, RV, RW, RO, TT;
// extern data UU, VV, WW;
// extern data FU, FV, FW, FR, FT;
// extern data ZRU, ZRV, ZRW, ZRO, ZTT;
// extern data WW1, WW2, WW3;
// extern data BX, BY, BZ;
// extern data ZBX, ZBY, ZBZ;
void print_DATA(data P) {
    for (int z = 0; z < nz; z++)
        for (int y = 0; y < ny; y++)
            for (int x = 0; x < nx; x++) {
                printi(x);
                printi(y);
                printi(z);
                printd(P[z][y][x]);
                endl();
            }
}
void print_DATAS() {
    fp = fopen("debug.datas", "w");
    print_DATA(RO);
    // print_DATA(RU); print_DATA(RV); print_DATA(RW); print_DATA(RO); print_DATA(TT);
    // print_DATA(UU); print_DATA(VV); print_DATA(WW);
    // print_DATA(FU); print_DATA(FV); print_DATA(FW); print_DATA(FR); print_DATA(FT);
    // print_DATA(ZRU); print_DATA(ZRV); print_DATA(ZRW); print_DATA(ZRO); print_DATA(ZTT);
    // print_DATA(WW1); print_DATA(WW2); print_DATA(WW3);
    // print_DATA(BX); print_DATA(BY); print_DATA(BZ);
    // print_DATA(ZBX); print_DATA(ZBY); print_DATA(ZBZ);
    fclose(fp);
}

/*AJACOBI*/
// extern double EXX[nx], dxxdx[nx], d2xxdx2[nx], DDX[nx];
// extern double WYY[ny], dyydy[ny], d2yydy2[ny], DDY[ny];
// extern double ZEE[nz], dzzdz[nz], d2zzdz2[nz], DDZ[nz];
void print_JACOBI() {
    // printf("XJacobi\n");
    for_i(nx) {
        printd(EXX[i]);
        printd(dxxdx[i]);
        printd(d2xxdx2[i]);
        printd(DDX[i]);
        endl();
    }
    // printf("YJacobi\n");
    for_i(ny) {
        printd(WYY[i]);
        printd(dyydy[i]);
        printd(d2yydy2[i]);
        printd(DDY[i]);
        endl();
    }
    // printf("ZJacobi\n");
    for_i(nz) {
        printd(ZEE[i]);
        printd(dzzdz[i]);
        printd(d2zzdz2[i]);
        printd(DDZ[i]);
        endl();
    }
}

/*ITER*/
// extern int nTotal, nStep0, nIt;
/*COMMUN*/
// extern int myPE, myPEy, myPEz;
// extern MPI_Datatype MPISIZE;
/*CPAR*/
// extern double CV, OCV, ORE, RE, REPR, THETA, GRAV, AMPT, GAMMA;
/*CROT*/
// extern double OMX, OMZ;
/*CMAG*/
// extern double ORM, RM, OBETA, AMPB, BFH, BZP;
/*CPEN*/
// extern double PZP, SIGMA, RKAPST, TB, RKAPM;
// extern double RKAPA[nz], DKAPA[nz];
void print_CPEN() {
    printd(PZP);
    printd(SIGMA);
    printd(RKAPST);
    printd(TB);
    printd(BZP);
    endl();
    for_i(nz) {
        printd(RKAPA[i]);
        printd(DKAPA[i]);
        endl();
    }
}
/*CPER*/
// extern double TP, XP, YP, ZP, TC, QFH, HH;
/*BOUNDS*/
// extern double xMax, yMax, zMax;
/*GRID*/
// extern double DD, hx, h2x, hy, h2y, hz, h2z, C13, C23, C43;
/*CTIM*/
// extern double DT, TIMT, TIMC, TIMI;
/*TRACE*/
// extern double UMACH;
/*RUNGKU*/
// extern double GAM1, GAM2, GAM3, ZETA1, ZETA2;
/*SPLINEX*/
// extern double HHX[nPx + 3], SIGX[nPx], AAX[nPx + 3], BBX[nPx + 3], XPINV[nPx], xHH;
// extern int iSEGx, kLOx[nPx + 3], kHIx[nPx + 3];
/*SPLINEY*/
// extern double HHY[nRy + 3], SIGY[nRy], AAY[nRy + 3], BBY[nRy + 3], YPINV[nRy], yHH;
// extern int iSEGy, kLOy[nRy + 3], kHIy[nRy + 3];
/*RELAX*/
// extern double TSTART, TOFF, RLAX;
/*SPECIALBOUND*/
// extern double TU, DZTB, DZTU;
/*FFCOM*/
// extern double A, C;


void print_with_tag_index(const data d, const char * tag, int Cx, int Cy, int Cz){
    /*
    print with tag and index, only first process will do this
    Cx: index in x with 0-index
    Cy: index in y with 0-index
    Cz: index in z with 0-index
    */
   if (myPE == 0){
        printf("%s[%d][%d][%d]=%21.14E\n", tag, Cz, Cy, Cx, d[Cz][Cy][Cx]);
   }
}

void print_with_tag_index_all(const char * tag, int Cx, int Cy, int Cz){
    /*
    print with tag and index, only first process will do this
    Cx: index in x with 0-index
    Cy: index in y with 0-index
    Cz: index in z with 0-index
    */
    /*
    if (myPE == 0){
        printf("current tag %s\n", tag);
        print_with_tag_index(RU, "RU", Cx, Cy, Cz);
        print_with_tag_index(RV, "RV", Cx, Cy, Cz);
        print_with_tag_index(RW, "RW", Cx, Cy, Cz);
        print_with_tag_index(RO, "RO", Cx, Cy, Cz);
        print_with_tag_index(TT, "TT", Cx, Cy, Cz);
        print_with_tag_index(UU, "UU", Cx, Cy, Cz);
        print_with_tag_index(VV, "VV", Cx, Cy, Cz);
        print_with_tag_index(WW, "WW", Cx, Cy, Cz);
        print_with_tag_index(FU, "FU", Cx, Cy, Cz);
        print_with_tag_index(FV, "FV", Cx, Cy, Cz);
        print_with_tag_index(FW, "FW", Cx, Cy, Cz);
        print_with_tag_index(FR, "FR", Cx, Cy, Cz);
        print_with_tag_index(FT, "FT", Cx, Cy, Cz);
        print_with_tag_index(ZRU, "ZRU", Cx, Cy, Cz);
        print_with_tag_index(ZRV, "ZRV", Cx, Cy, Cz);
        print_with_tag_index(ZRW, "ZRW", Cx, Cy, Cz);
        print_with_tag_index(ZRO, "ZRO", Cx, Cy, Cz);
        print_with_tag_index(ZTT, "ZTT", Cx, Cy, Cz);
        print_with_tag_index(WW1, "WW1", Cx, Cy, Cz);
        print_with_tag_index(WW2, "WW2", Cx, Cy, Cz);
        print_with_tag_index(WW3, "WW3", Cx, Cy, Cz);
        print_with_tag_index(BX, "BX", Cx, Cy, Cz);
        print_with_tag_index(BY, "BY", Cx, Cy, Cz);
        print_with_tag_index(BZ, "BZ", Cx, Cy, Cz);
        print_with_tag_index(ZBX, "ZBX", Cx, Cy, Cz);
        print_with_tag_index(ZBY, "ZBY", Cx, Cy, Cz);
        print_with_tag_index(ZBZ, "ZBZ", Cx, Cy, Cz);
   }
   */
}

#endif
