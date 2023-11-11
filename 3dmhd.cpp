/**
 * This is a finite difference program for the Cray T3E to simulate three
 * dimensional magnetohydrodynamics in a rotating fully ionized
 * hydrogen ideal-gas of adiabatic stratification in the upper potion and
 * stable stratification in the lower.  Rotation and magnetic fields can be
 * turned on and off with logical parameters lRot and lMag.  The code empolys
 * irregular grids, second-order finite differencing, and fully explicit
 * time stepping.  Must be linked with 3dmhdset.f and 3dmhdsub.f.
 * Unless noted the code was written by Mark Peter Rast.
 * Original two-dimensional version finished 12/27/95.  MPP version 11/1/96.
 * Three-dimensional magnetohydrodynamic version 2/5/98.
 * Two dimensional decomposition of the domain by Matthias Rempel 03/21/02.
 */
#include "3dmhdio.h"
#include "3dmhdparam.h"
#include "3dmhdphsic.h"
#include "3dmhdprint.h"
#include "3dmhdset.h"
#include "3dmhdsub.h"
#include "3dmhdutils.h"

#include <mpi.h>

typedef double data[nz][ny][nx];

int ipar[32];
double par[64];
/*BIG*/
data RU, RV, RW, RO, TT;
data UU, VV, WW;
data FU, FV, FW, FR, FT;
data ZRU, ZRV, ZRW, ZRO, ZTT;
data WW1, WW2, WW3;
data BX, BY, BZ;
data ZBX, ZBY, ZBZ;
/*AJACOBI*/
double EXX[nx], dxxdx[nx], d2xxdx2[nx], DDX[nx];
double WYY[ny], dyydy[ny], d2yydy2[ny], DDY[ny];
double ZEE[nz], dzzdz[nz], d2zzdz2[nz], DDZ[nz];
/*ITER*/
int nTotal, nStep0, nIt;
/*COMMUN*/
int myPE, myPEy, myPEz;
MPI_Datatype MPISIZE;
/*CPAR*/
double CV, OCV, ORE, RE, REPR, THETA, GRAV, AMPT, GAMMA;
/*CROT*/
double OMX, OMZ;
/*CMAG*/
double ORM, RM, OBETA, AMPB, BFH, BZP;
/*CPEN*/
double PZP, SIGMA, RKAPST, TB, RKAPM;
double RKAPA[nz], DKAPA[nz];
/*CPER*/
double TP, XP, YP, ZP, TC, QFH, HH;
/*BOUNDS*/
double xMax, yMax, zMax;
/*GRID*/
double DD, hx, h2x, hy, h2y, hz, h2z;

double C13, C23, C43;
/*CTIM*/
double DT, TIMT, TIMC, TIMI;
/*TRACE*/
double UMACH;
/*RUNGKU*/
double GAM1, GAM2, GAM3, ZETA1, ZETA2;
/*SPLINEX*/
double HHX[nPx + 3], SIGX[nPx], AAX[nPx + 3], BBX[nPx + 3], XPINV[nPx], xHH;
int kLOx[nPx + 3], kHIx[nPx + 3], iSEGx;
/*SPLINEY*/
double HHY[nRy + 3], SIGY[nRy], AAY[nRy + 3], BBY[nRy + 3], YPINV[nRy], yHH;
int kLOy[nRy + 3], kHIy[nRy + 3], iSEGy;
/*RELAX*/
double TSTART, TOFF, RLAX;
/*SPECIALBOUND*/
double TU, DZTB, DZTU;
/*FFCOM*/
double A, C;

void init_spline_x(double* EX);
void init_spline_y(double* EY);
std::string PEextn();
void set();

int main(int argc, char* argv[]) {
    int NN;
    auto rt = mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (!rt) printf("Folder created\n");else printf("Impossible create folder ");
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE);
    MPI_Comm_size(MPI_COMM_WORLD, &NN);
    if (NN != nPE) exit_err("MAIN:  Number of processors not equal to nPE");

    myPEy = myPE % nPEy;
    myPEz = myPE / nPEy;
    printf("[%02d] myPEy [%02d] myPEz [%02d]\n", myPE, myPEy, myPEz);

    if (iWord == 8) {
        MPISIZE = MPI_DOUBLE_PRECISION;
    } else {
        MPISIZE = MPI_FLOAT;
    }

    set();

    // TODO 以上的部分验证完毕

    std::string fIn0 = fInp + ".dat0." + PEextn();
    std::string fIn1 = fInp + ".par";

    if (nStart == 0) {
        Static();
        TIMI = 0.0f;
    } else {
#if lRem
        Static();
#endif

        if (myPE == 0) loadParam(fIn1, 0);
        MPI_Bcast(par1, 96, MPISIZE, 0, MPI_COMM_WORLD);
        int nx1 = static_cast<int>(par1[10]);
        int ny1 = static_cast<int>(par1[11]) / static_cast<int>(par1[20]);
        int nz1 = static_cast<int>(par1[12]) / static_cast<int>(par1[13] / par1[20]);
        if (nx - ix != nx1 || ny - iy != ny1 || nz != nz1 + iz) exit_err("MAIN:  Grid size mismatch, no interpolation");

        if (myPE == 0) loadParam(fIn1, nStart - 1);
        MPI_Bcast(par1, 96, MPISIZE, 0, MPI_COMM_WORLD);
        double TIMI = par1[63] + par1[64];
        loadData(fIn0, nStart - 1);
        if (myPEz == nPEz - 1) {
            TB = TT[1][1][nz - iz / 2 - 1];
        }
    }

    // MPI_Barrier(MPI_COMM_WORLD);
    if (myPE == 0) print_DATAS();
    // MPI_Barrier(MPI_COMM_WORLD);
    // exit_err("debug");
    // TODO check

    communicate();

    CV = 1.0f / (GAMMA - 1.0f);
    OCV = 1.0f / CV;
    for (int z = 0; z < nz; z++) {
        for (int y = 1; y < ny - iy + 1; y++) {
            for (int x = 1; x < nx - ix + 1; x++) {
                UU[z][y][x] = (square(RU[z][y][x]) + square(RV[z][y][x]) + square(RW[z][y][x])) / square(RO[z][y][x]);
                VV[z][y][x] = GAMMA * TT[z][y][x];
#if lMag
                VV[z][y][x] += OBETA * (square(BX[z][y][x]) + square(BY[z][y][x]) + square(BZ[z][y][x])) / square(RO[z][y][x]);
#endif
                UU[z][y][x] += VV[z][y][x] + 2.0f * sqrt(UU[z][y][x] * VV[z][y][x]);
            }
        }
    }
#define ISW 1
#if ISW == 0
    double RMIN = minVal(RO, 0, nz, 1, ny - iy + 1, 1, nx - ix + 1);
    double VMAX = sqrt(maxVal(UU, 0, nz, 1, ny - iy + 1, 1, nx - ix + 1));

    if (nx > ix + 1) {
        double minDDX = minVal(DDX, 1, nx - ix + 1);
        double minDDY = minVal(DDY, 1, ny - iy + 1);
        double minDDZ = minVal(DDZ, 0, nz);
        DD = min(min(minDDX, minDDY), minDDZ);
    } else {
        double minDDY = minVal(DDY, 1, ny - iy + 1);
        double minDDZ = minVal(DDZ, 0, nz);
        DD = min(minDDY, minDDZ);
    }

    RKAPM = maxVal(RKAPA, 0, nz);

    wMin[0] = RMIN;
    wMin[1] = DD;
    MPI_Allreduce(wMin, wMout, 2, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    RMIN = wMout[0];
    DD = wMout[1];

    wMin[0] = VMAX;
    wMin[1] = RKAPM;
    MPI_Allreduce(wMin, wMout, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    VMAX = wMout[0];
    RKAPM = wMout[1];

    DT1 = DD / VMAX;

    #if lRem
    DT2 = 0.5f * DD * DD * REPR * CV * RMIN / (1.0f + RKAPM);
    #else
    DT2 = 0.5f * DD * DD * REPR * CV * RMIN / RKAPM;
    #endif

    DT3 = 0.375f * DD * DD * RE * RMIN;

    #if lMag
    DT4 = 0.5f * DD * DD * RM;
    DT = SFF * min(min(DT1, DT2), min(DT3, DT4));
    #else
    DT4 = 0.0f;
    DT = SFF * min(min(DT1, DT2), DT3);
    #endif
#else
    for (int z = iz / 2; z < nz - iz / 2; z++) {
        for (int y = 1; y < ny - iy + 1; y++) {
            for (int x = 1; x < nx - ix + 1; x++) {
                if (nx > ix + 1) {
                    DD = min(min(DDX[x], DDY[y]), DDZ[z]);
                } else {
                    DD = min(DDY[y], DDZ[z]);
                }
                WW1[z][y][x] = DD / sqrt(UU[z][y][x]);
    #if lRem
                WW2[z][y][x] = 0.5f * DD * DD * REPR * CV * RO[z][y][x] / (1.0f + RKAPA[z]);
    #else
                WW2[z][y][x] = 0.5f * DD * DD * REPR * CV * RO[z][y][x] / RKAPA[z];
    #endif
                WW3[z][y][x] = 0.375f * DD * DD * RE * RO[z][y][x];
    #if lMag
                VV[z][y][x] = 0.5f * DD * DD * RM;
    #endif
            }
        }
    }

    wMin[0] = minVal(WW1, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
    wMin[1] = minVal(WW2, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
    wMin[2] = minVal(WW3, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);

    int MINCNT = 3;

    #if lMag
    wMin[3] = minVal(VV, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
    MINCNT = 4;
    #endif

    MPI_Allreduce(wMin, wMout, MINCNT, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    DT1 = wMout[0];
    DT2 = wMout[1];
    DT3 = wMout[2];

    #if lMag
    DT4 = wMout[3];
    DT = SFF * min(min(DT1, DT2), min(DT3, DT4));
    #else
    DT4 = 0.0f;
    DT = SFF * min(min(DT1, DT2), DT3);
    #endif
#endif
#undef ISW

    fill(FU, 0, 0, nz, 0, ny, 0, nx);
    fill(FV, 0, 0, nz, 0, ny, 0, nx);
    fill(FW, 0, 0, nz, 0, ny, 0, nx);
    fill(FT, 0, 0, nz, 0, ny, 0, nx);
    fill(FR, 0, 0, nz, 0, ny, 0, nx);
    fill(WW1, 0, 0, nz, 0, ny, 0, nx);
    fill(WW2, 0, 0, nz, 0, ny, 0, nx);
    fill(WW3, 0, 0, nz, 0, ny, 0, nx);

    std::string ff0 = fOut + ".dat0." + PEextn();
    std::string ff1 = fOut + ".par";
    std::string fLog = fOut + ".lis";

    for (int i = 0; i < 80; i++) blanks[i] = ' ';
    for (int i = 0; i < 32; i++) par1[i] = float(ipar[i]);
    for (int i = 0; i < 64; i++) par1[i + 32] = par[i];

    if (nIt == 0 && myPE == 0) logGreet(fLog);

    if (myPE == 0) printf("iteration  0: DT = %20.14E\n", DT);

    int nBeg = nIt + 1;
    for (int nk = nBeg; nk <= nTotal; nk++) {
        step();
        par1[14] = nIt;
        par1[62] = DT;
        if (nk == nTotal) {
            par1[63] = TIMT;
            par1[64] = 0.0f;
        } else {
            par1[63] = TIMI;
            par1[64] = TIMC;
        }
        if (myPE == 0) printf("iteration %2d: DT = %20.14E\n", nk, DT);

        if (DT < 1.0e-8f) {
            dump(ff0, ff1);

            for (int z = iz / 2; z < nz - iz / 2; z++) {
                for (int y = 1; y < ny - iy + 1; y++) {
                    for (int x = 1; x < nx - ix + 1; x++) {
                        DD = min(min(DDX[x], DDY[y]), DDZ[z]);
                        WW1[z][y][x] = DD / sqrt(UU[z][y][x]);
#if lRem
                        WW2[z][y][x] = 0.5f * DD * DD * REPR * CV * RO[z][y][x] / (1.0f + RKAPA[z]);
#else
                        WW2[z][y][x] = 0.5f * DD * DD * REPR * CV * RO[z][y][x] / RKAPA[z];
#endif
                        WW3[z][y][x] = 0.375f * DD * DD * RE * RO[z][y][x];
#if lMag
                        VV[z][y][x] = 0.5f * DD * DD * RM;
#endif
                    }
                }
            }
            double min_val;
            int min_pos[3];
            min_val = minVal(WW1, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
            minLoc(WW1, min_pos, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
            printf("[%02d] Min DT1: %10.4E %d %d %d\n", myPE, min_val, min_pos[0], min_pos[1], min_pos[2]);
            min_val = minVal(WW2, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
            minLoc(WW2, min_pos, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
            printf("[%02d] Min DT2: %10.4E %d %d %d\n", myPE, min_val, min_pos[0], min_pos[1], min_pos[2]);
            min_val = minVal(WW3, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
            minLoc(WW3, min_pos, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
            printf("[%02d] Min DT3: %10.4E %d %d %d\n", myPE, min_val, min_pos[0], min_pos[1], min_pos[2]);
#if lMag
            min_val = minVal(VV, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
            minLoc(VV, min_pos, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
            printf("[%02d] Min DT4: %10.4E %d %d %d\n", myPE, min_val, min_pos[0], min_pos[1], min_pos[2]);
#endif

            if (myPE == 0) logFatal(fLog);
            exit_err("MAIN  DT too small, halt");
        }

        if (nIt % nStep0 == 0) {
            dump(ff0, ff1);
            if (myPE == 0) printf("iteration %5d of %5d\n", nIt, nTotal);
            if (myPE == 0) logEpoch(fLog);
        }
    }

    if (myPE == 0) {
        std::ofstream outLog(fLog, std::ios::out | std::ios::app);
        s4000(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
        s4090(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
        s4000(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
        s4000(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
        outLog.close();
    }

    MPI_Finalize();
    return 0;
}

void set() {
    C13 = 1.0f / 3.0f;
    C23 = 2.0f * C13;
    C43 = 4.0f * C13;

    GAM1 = 8.0f / 15.0f;
    GAM2 = 5.0f / 12.0f;
    GAM3 = 3.0f / 4.0f;
    ZETA1 = -17.0f / 60.0f;
    ZETA2 = -5.0f / 12.0f;

    hx = 0.5f * float(nPx - 1);
    h2x = float(nPx - 1) * float(nPx - 1);
    hy = 0.5f * float(nPy - 1);
    h2y = float(nPy - 1) * float(nPy - 1);
    hz = 0.5f * float(nPz - 1);
    h2z = float(nPz - 1) * float(nPz - 1);

    setup(ipar, par);

    nCase = ipar[0];
    nCaseP = ipar[1];
    nTotal = ipar[2];
    nStep0 = ipar[3];
    nStart = ipar[4];

    nIt = 0;
    nDump0 = 0;
    NSW = 0;

    ipar[5] = ixCon;
    ipar[6] = iyCon;
    ipar[7] = izCon;
    ipar[8] = iTCon;
    ipar[9] = iBCon;

    ipar[10] = nPx;
    ipar[11] = nPy;
    ipar[12] = nPz;
    ipar[13] = nPE;

    ipar[15] = nTube;
    ipar[16] = nGrid;
    ipar[17] = nCol;
    ipar[18] = nRow;
    ipar[19] = iD;
    ipar[20] = nPEy;

    RE = par[0];
    PR = par[1];
    THETA = par[2];
    GRAV = par[3];
    RY = par[4];
    ANG = par[5];
    RM = par[6];
    BETA = par[7];
    PZP = par[8];
    SIGMA = par[9];
    POLYS = par[10];
    TP = par[11];
    TC = par[12];
    XP = par[13];
    YP = par[14];
    ZP = par[15];
    xMax = par[16];
    yMax = par[17];
    zMax = par[18];
    AMPT = par[19];
    AMPB = par[20];
    BFH = par[21];
    BZP = par[22];
    QFH = par[52];
    GAMMA = par[53];
    HH = par[54];

    TSTART = par[57];
    TOFF = par[58];
    RLAX = par[59];

#if lRot
    OMX = sin(ANG) / RY;
    OMZ = cos(ANG) / RY;
#else
    OMX = 0.0f;
    OMZ = 0.0f;
#endif

    par[23] = XX1, par[24] = XX2;
    par[25] = YY1, par[26] = YY2;
    par[27] = ZZ1, par[28] = ZZ2;
    par[29] = SFF;

    par[38] = XA, par[39] = XB, par[40] = XC, par[41] = XD;
    par[42] = YA, par[43] = YB, par[44] = YC, par[45] = YD;
    par[46] = ZA, par[47] = ZB, par[48] = ZC, par[49] = ZD;

    par[55] = DH, par[56] = DV;

    REPR = RE * PR;
    ORE = 1.0f / RE;
    if (GRAV == 0.0) {
        RKAPST = 1.0f;
    } else {
        RKAPST = (POLYS + 1.0f) * THETA / GRAV;
    }
    ORM = 1.0f / RM;
    OBETA = 2.0f / BETA;
    TIMC = 0.0f;
    if (nx > ix + 1) {
        xJacobi();
    } else {
        fill(EXX, 0, nx);
        fill(dxxdx, 0, nx);
        fill(d2xxdx2, 0, nx);
        fill(DDX, 0, nx);
    }
    yJacobi();
    zJacobi();

    if (nx > ix + 1) {
        init_spline_x(&EXX[ix - 1]);
    }
    init_spline_y(&WYY[iy - 1]);

#if !lRem
    if (PZP == 0.0) {
        for (int i = 0; i < nz; i++) {
            RKAPA[i] = 1.0f;
            DKAPA[i] = 0;
        }
    } else {
        for (int i = 0; i < nz; i++) {
            RKAPA[i] = 1.0f + (RKAPST - 1.0f) / 2.0f * (1.0f + tanh((ZEE[i] - PZP) / SIGMA));
            DKAPA[i] = (RKAPST - 1.0f) / 2.0f / SIGMA / square(cosh((ZEE[i] - PZP) / SIGMA));
        }
    }
#endif

    // MPI_Barrier(MPI_COMM_WORLD);
    // if (myPE == 0) print_CPEN();
    // MPI_Barrier(MPI_COMM_WORLD);
    // exit_err("debug");
}

std::string PEextn() {
    std::string str(4, '0');
    for (int k = 3, t = myPE; k >= 0; k--)
        str[k] = '0' + t % 10, t /= 10;
    return str;
}

void init_spline_x(double* EX) {
    double X2[nPx];
    double EX2[nPx + 3];
    iSEGx = nPx - 1;
    while ((iSEGx % 4) != 0) iSEGx += 1;
    double xMin = minVal(EX, 0, nPx);
    double xMax = maxVal(EX, 0, nPx);
    xHH = (xMax - xMin) / float(iSEGx);
    for (int i = 0; i <= iSEGx; i++) EX2[i] = xHH * float(i) + xMin;
    X2[0] = 0.0f;
    for (int i = 1; i < nPx - 1; i++) {
        SIGX[i] = (EX[i] - EX[i - 1]) / (EX[i + 1] - EX[i - 1]);
        XPINV[i] = 1.0f / (SIGX[i] * X2[i - 1] + 2.0f);
        X2[i] = (SIGX[i] - 1.0f) * XPINV[i];
    }
    for (int i = 0; i <= iSEGx; i++) {
        double X = EX2[i];
        kLOx[i] = 0;
        kHIx[i] = nPx - 1;
        while (kHIx[i] - kLOx[i] > 1) {
            int k = (kHIx[i] + kLOx[i]) / 2;
            if (EX[k] > X) {
                kHIx[i] = k;
            } else {
                kLOx[i] = k;
            }
        }
        HHX[i] = EX[kHIx[i]] - EX[kLOx[i]];
        if (HHX[i] == 0) exit_err("init_spline_x: EX values must be distinct.");
        AAX[i] = (EX[kHIx[i]] - X) / HHX[i];
        BBX[i] = (X - EX[kLOx[i]]) / HHX[i];
    }
}

void init_spline_y(double* EY) {
    double Y2[nRy];
    double EY2[nRy + 3];
    iSEGy = nRy - 1;
    while ((iSEGy % 4) != 0) iSEGy += 1;
    double yMin = minVal(EY, 0, nRy);
    double yMax = maxVal(EY, 0, nRy);
    yHH = (yMax - yMin) / float(iSEGy);
    for (int i = 0; i <= iSEGy; i++) EY2[i] = yHH * float(i) + yMin;
    Y2[0] = 0.0f;
    for (int i = 1; i < nRy - 1; i++) {
        SIGY[i] = (EY[i] - EY[i - 1]) / (EY[i + 1] - EY[i - 1]);
        XPINV[i] = 1.0f / (SIGY[i] * Y2[i - 1] + 2.0f);
        Y2[i] = (SIGY[i] - 1.0f) * XPINV[i];
    }
    for (int i = 0; i <= iSEGy; i++) {
        double Y = EY2[i];
        kLOy[i] = 0;
        kHIy[i] = nRy - 1;
        while (kHIy[i] - kLOy[i] > 1) {
            int k = (kHIy[i] + kLOy[i]) / 2;
            if (EY[k] > Y) {
                kHIy[i] = k;
            } else {
                kLOy[i] = k;
            }
        }
        HHY[i] = EY[kHIy[i]] - EY[kLOy[i]];
        if (HHY[i] == 0) exit_err("init_spline_y: EY values must be distinct.");
        AAY[i] = (EY[kHIy[i]] - Y) / HHY[i];
        BBY[i] = (Y - EY[kLOy[i]]) / HHY[i];
    }
}
