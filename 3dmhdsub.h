#ifndef _3DMHD_SUB_H_
#define _3DMHD_SUB_H_

#include "3dmhddata.h"
#include "3dmhdmpi.h"
#include "3dmhdparam.h"
#include "3dmhdutils.h"
#define VIEW_INDX 10,10,10

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>

void fluxes(FORRANGEPARAM);
void BCON(FORRANGEPARAM);
void step_block_GAM(double COEF, FORRANGEPARAM);
void step_block_ZETA(double COEF, FORRANGEPARAM);
void step_block_ZETA0(FORRANGEPARAM);
void horizontal_mean(double varM[nz], data var);

void dump(data var, std::string filename) {
    std::ofstream file("./dump/" + filename + ".bin", std::ios::binary);
    file.write((char*)var, sizeof(double) * nx * ny * nz);
    file.close();
}

void step() {
    double WMin[4], WMout[4];
    const int zl = iz / 2, zr = nz - iz / 2;
    const int yl = iy / 2, yr = ny - iy / 2;
    const int xl = ix / 2, xr = nx - ix / 2;
    UMACH = 0.0f;
    double RMIN = 1e9f;
    double VMAX = 0.0f;
    double OGAMMA = 1.0f / GAMMA;
    int z, y, x;
    FOR3D UU[z][y][x] = (RU[z][y][x] * RU[z][y][x] + RV[z][y][x] * RV[z][y][x] + RW[z][y][x] * RW[z][y][x]) / (RO[z][y][x] * RO[z][y][x]);
    UMACH = 0.0f;
    FOR3D UMACH = max(UMACH, OGAMMA * UU[z][y][x] / TT[z][y][x]);
    UMACH = sqrt(UMACH);
#if lMag
    FOR3D UU[z][y][x] += GAMMA * TT[z][y][x] + (BX[z][y][x] * BX[z][y][x] + BY[z][y][x] * BY[z][y][x] + BZ[z][y][x] * BZ[z][y][x]) * OBETA / RO[z][y][x] + 2.0f * sqrt(UU[z][y][x] * (GAMMA * TT[z][y][x] + (BX[z][y][x] * BX[z][y][x] + BY[z][y][x] * BY[z][y][x] + BZ[z][y][x] * BZ[z][y][x]) * OBETA / RO[z][y][x]));
#else
    FOR3D UU[z][y][x] += GAMMA * TT[z][y][x] + 2.0f * sqrt(GAMMA * UU[z][y][x] * TT[z][y][x]);
#endif
    FOR3D VMAX = max(VMAX, UU[z][y][x]), RMIN = min(RMIN, RO[z][y][x]);
    VMAX = sqrt(VMAX);
    MPI_Allreduce(&UMACH, WMout, 1, MPISIZE, MPI_MAX, MPI_COMM_WORLD);
    UMACH = WMout[0];

#define ISW 1
#if ISW == 0
    WMin[0] = RMIN;
    WMin[1] = -VMAX;
    MPI_Allreduce(WMin, WMout, 2, MPISIZE, MPI_MIN, MPI_COMM_WORLD);
    RMIN = WMout[0];
    VMAX = -WMout[1];
    double DT1 = DD / VMAX, DT2;
    #if lRem
    DT2 = 0.5f * DD * DD * REPR * CV * RMIN / (1.0f + RKAPM);
    #else
    DT2 = 0.5f * DD * DD * REPR * CV * RMIN / RKAPM;
    #endif
    double DT3 = 0.375f * DD * DD * RE * RMIN;
    #if lMag
    DT = SFF * min(min(DT1, DT2), min(DT3, 0.5f * DD * DD * RM));
    #else
    DT = SFF * min(min(DT1, DT2), DT3);
    #endif
#else
    FOR3D {
        DD = min(DDY[y], DDZ[z]);
        if (nx > ix + 1) DD = min(DD, DDX[x]);
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
    WMin[0] = minVal(WW1, FORRANGE);
    WMin[1] = minVal(WW2, FORRANGE);
    WMin[2] = minVal(WW3, FORRANGE);
    // above correct
    // if (myPE == 0) {
    //     printf("WMin[0]=%20.14E\nWMin[1]=%20.14E\nWMin[2]=%20.14E\n", WMin[0], WMin[1], WMin[2]);
    // }
    int MINCNT = 3;
    #if lMag
    WMin[3] = minVal(VV, FORRANGE), MINCNT = 4;
    #endif
    MPI_Allreduce(WMin, WMout, MINCNT, MPISIZE, MPI_MIN, MPI_COMM_WORLD);
    #if lMag
    DT = SFF * min(min(WMout[0], WMout[1]), min(WMout[2], WMout[3]));
    #else
    DT = SFF * min(min(WMout[0], WMout[1]), WMout[2]);
    #endif
#endif
#undef ISW

    /**
     * Third-order Runge-Kutta timestepping scheme of Wray (Spalart, P.R.,
     * Moser, R.D., & Rogers M.M., Spectral Methods for the Navier-Stokes
     * Equations with One Infinite and Two Periodic Directions,
     * j. Comp. Phys., 96 297-324 1991).
     **/
    // Calculate the first Runge-Kutta substep.
    step_block_ZETA0(FORRANGE);
    fluxes(FORRANGE);
    print_with_tag_index_all("after fluxes1", VIEW_INDX);
    step_block_GAM(GAM1 * DT, FORRANGE);
    print_with_tag_index_all("after fluxes1 GAM", VIEW_INDX);

    BCON(FORRANGE);
    print_with_tag_index_all("after BCON1", VIEW_INDX);
    communicate();
    print_with_tag_index_all("after communicate1", VIEW_INDX);
    // Calculate second Runge-Kutta substep.
    step_block_ZETA(ZETA1 * DT, FORRANGE);
    print_with_tag_index_all("step_block_ZETA 1", VIEW_INDX);
    fluxes(FORRANGE);
    print_with_tag_index_all("after fluxes2", VIEW_INDX);
    step_block_GAM(GAM2 * DT, FORRANGE);
    print_with_tag_index_all("after fluxes2 GAM", VIEW_INDX);
    BCON(FORRANGE);
    print_with_tag_index_all("after BCON2", VIEW_INDX);
    communicate();
    print_with_tag_index_all("after communicate2", VIEW_INDX);
    // Calculate third Runge-Kutta substep.
    step_block_ZETA(ZETA2 * DT, FORRANGE);
    print_with_tag_index_all("step_block_ZETA 2", VIEW_INDX);
    fluxes(FORRANGE);
    print_with_tag_index_all("after fluxes3", VIEW_INDX);
    step_block_GAM(GAM3 * DT, FORRANGE);
    print_with_tag_index_all("after fluxes3 GAM", VIEW_INDX);
    BCON(FORRANGE);
    print_with_tag_index_all("after BCON3", VIEW_INDX);
    communicate();
    print_with_tag_index_all("after communicate3", VIEW_INDX);


    nIt += 1;
    TIMC += DT;
    TIMT = TIMI + TIMC;
}

void step_block_ZETA0(FORRANGEPARAM) {
    copy(ZRU, RU, FORRANGE);
    copy(ZRV, RV, FORRANGE);
    copy(ZRW, RW, FORRANGE);
    copy(ZRO, RO, FORRANGE);
    copy(ZTT, TT, FORRANGE);
#if lMag
    copy(ZBX, BX, FORRANGE);
    copy(ZBY, BY, FORRANGE);
    copy(ZBZ, BZ, FORRANGE);
#endif
}

void step_block_ZETA(double COEF, FORRANGEPARAM) {
    scale_add(ZRU, RU, FU, COEF, FORRANGE);
    scale_add(ZRV, RV, FV, COEF, FORRANGE);
    scale_add(ZRW, RW, FW, COEF, FORRANGE);
    scale_add(ZRO, RO, FR, COEF, FORRANGE);
    scale_add(ZTT, TT, FT, COEF, FORRANGE);
#if lMag
    scale_add(ZBX, BX, WW1, COEF, FORRANGE);
    scale_add(ZBY, BY, WW2, COEF, FORRANGE);
    scale_add(ZBZ, BZ, WW3, COEF, FORRANGE);
#endif
}

void step_block_GAM(double COEF, FORRANGEPARAM) {
    scale_add(RU, ZRU, FU, COEF, FORRANGE);
    scale_add(RV, ZRV, FV, COEF, FORRANGE);
    scale_add(RW, ZRW, FW, COEF, FORRANGE);
    scale_add(RO, ZRO, FR, COEF, FORRANGE);
    scale_add(TT, ZTT, FT, COEF, FORRANGE);
#if lMag
    scale_add(BX, ZBX, WW1, COEF, FORRANGE);
    scale_add(BY, ZBY, WW2, COEF, FORRANGE);
    scale_add(BZ, ZBZ, WW3, COEF, FORRANGE);
#endif
}

void fluxes(FORRANGEPARAM) {
    static int ICALL = 0;
    static double TTM_MIN = 1e99;
    static double TTM_MAX = -1e99;
    static double TFAC = 1.0;

    static double TTM[nz], ROM[nz], HRAD[nz], FCONM[nz], FRADM[nz], ALPHA[nz], CORRECT[nz], ADDSUM[nPE], ENDVAL[nPE], WWY[ny], wwz[nz];

    double tmp;

#if ixCon == 0 && iyCon == 0
    func_d1<OPT::SET, DIM::X>(FORRANGE, FR, RU, -1); // do 10.1
    func_d1<OPT::SUB, DIM::Y>(FORRANGE, FR, RV);     // do 10.2
#else
    exit_err("fluxes: Non-periodic horizontal boundaries");
#endif
    func_d1<OPT::SET, DIM::Z>(FORRANGE, WW1, RW); // do 40
    if (myPEz == 0) {
        double TMPZ = hz * dzzdz[zl];
        for (int y = yl; y < yr; y++) {
            for (int x = xl; x < xr; x++) {
                WW1[zl][y][x] = (4.0f * RW[zl + 1][y][x] - RW[zl + 2][y][x]) * TMPZ; // do 50
            }
        }
    }
    if (myPEz == nPEz - 1) {
        double TMPZ = hz * dzzdz[zr - 1];
        for (int y = yl; y < yr; y++) {
            for (int x = xl; x < xr; x++) {
                WW1[zr - 1][y][x] = (RW[zr - 3][y][x] - 4.0f * RW[zr - 2][y][x]) * TMPZ; // do 60
            }
        }
    }
    func_base<OPT::SUB>(FORRANGE, FR, WW1); // do 70
    print_with_tag_index_all("after fluxes do 70", VIEW_INDX);
#if iD != 0 || lRem
    if (myPEz == 0) {
        for (int z = 0; z < iz / 2; z++) {
            memcpy(&RO[z][0][0], &RO[z + iz / 2 + 1][0][0], sizeof(double) * nx * ny);
        }
    }
    if (myPEz == nPEz - 1) {
        for (int z = nz - iz / 2; z < nz; z++) {
            memcpy(&RO[z][0][0], &RO[z - iz / 2 - 1][0][0], sizeof(double) * nx * ny);
        }
    }
    horizontal_mean(ROM, RO);
    func_base<OPT::SET>(FORRANGE_EXT, WW1, RO);                 // do 9120.1
    func_base_extend<OPT::SUB, DIM::Z>(FORRANGE_EXT, WW1, ROM); // do 9120.2
#endif
#if iD != 0
    #if iD == 1
    for (int z = 0; z < nz; z++) wwz[z] = 1.0f / ROM[z];
    #endif
    #if iD == 2
    double SGM = 1.0f;
    double CLN = -4.0f * log(2.0f) / SGM / SGM;
    for (int z = 0; z < nz; z++) {
        wwz[z] = exp(CLN * square(ZEE[z])) + exp(CLN * square(ZEE[z] - zMax));
    }
    #endif
    #if DH != 0.0
    func_d2<OPT::SET, DIM::X>(FORRANGE, WW2, WW1);                       // do 9130
    func_d2<OPT::SET, DIM::Y>(FORRANGE, WW3, WW1);                       // do 9135
    func_mul_extend<OPT::ADD, DIM::Z>(FORRANGE, FR, WW2, wwz, DH * ORE); // do 9140.1
    func_mul_extend<OPT::ADD, DIM::Z>(FORRANGE, FR, WW3, wwz, DH * ORE); // do 9140.2
    #endif
    #if DV != 0.0
    func_d2<OPT::SET, DIM::Z>(FORRANGE, WW2, WW1);                       // do 9145
    func_mul_extend<OPT::ADD, DIM::Z>(FORRANGE, FR, WW2, wwz, DV * ORE); // do 9150
    #endif
#endif
    func_base<OPT::SET>(FORRANGE, FW, RO, GRAV);    // do 80
    func_mul2<OPT::SET>(FORRANGE_EXT, WW1, RO, TT); // do 90.1
    func_inv(FORRANGE_EXT, RO);                     // do 90.2

    func_d1<OPT::SUB, DIM::Z>(FORRANGE, FW, WW1);     // do 120.1
    func_d1<OPT::SET, DIM::Y>(FORRANGE, FV, WW1, -1); // do 120.2
    func_d1<OPT::SET, DIM::X>(FORRANGE, FU, WW1, -1); // do 120.3

    func_mul2<OPT::SET>(FORRANGE_EXT, UU, RU, RO); // do 160
    func_mul2<OPT::SET>(FORRANGE_EXT, VV, RV, RO); // do 170
    func_mul2<OPT::SET>(FORRANGE_EXT, WW, RW, RO); // do 180
    print_with_tag_index_all("after fluxes do 180", VIEW_INDX);

    func_mul_d1<OPT::SUB, DIM::X>(FORRANGE, FU, RU, UU); // do 210.1
    func_mul_d1<OPT::SUB, DIM::Y>(FORRANGE, FU, RU, VV); // do 210.2
    func_mul_d1<OPT::SUB, DIM::Z>(FORRANGE, FU, RU, WW); // do 210.3

    func_mul_d1<OPT::SUB, DIM::X>(FORRANGE, FV, RV, UU); // do 280.1
    func_mul_d1<OPT::SUB, DIM::Y>(FORRANGE, FV, RV, VV); // do 280.2
    func_mul_d1<OPT::SUB, DIM::Z>(FORRANGE, FV, RV, WW); // do 280.3

    func_mul_d1<OPT::SUB, DIM::X>(FORRANGE, FW, RW, UU); // do 330.1
    func_mul_d1<OPT::SUB, DIM::Y>(FORRANGE, FW, RW, VV); // do 330.2
    func_mul_d1<OPT::SUB, DIM::Z>(FORRANGE, FW, RW, WW); // do 330.3
    print_with_tag_index_all("after fluxes do 330", VIEW_INDX);


#if !lShr
    func_d1_mul<OPT::SET, DIM::X>(FORRANGE, FT, TT, UU, -1); // do 450.1
    func_d1_mul<OPT::SUB, DIM::Y>(FORRANGE, FT, TT, VV);     // do 450.2
    func_d1_mul<OPT::SUB, DIM::Z>(FORRANGE, FT, TT, WW);     // do 450.3
    if ((RLAX > 0.0) || lRem) {
        horizontal_mean(TTM, TT);
        int z, y, x;
        For(z, 0, nz) For(y, 0, ny) For(x, 0, nx) WW1[z][y][x] = TT[z][y][x] - TTM[z]; // do anonymous
    }
    print_with_tag_index_all("horizontal_mean", VIEW_INDX);
#endif
#if lRem
    tmp = OCV / REPR;
    func_d2_mul<OPT::ADD, DIM::X>(FORRANGE, FT, WW1, RO, tmp); // do 518.1
    func_d2_mul<OPT::ADD, DIM::Y>(FORRANGE, FT, WW1, RO, tmp); // do 518.2
    func_d2_mul<OPT::ADD, DIM::Z>(FORRANGE, FT, WW1, RO, tmp); // do 518.3
    for (int z = zl; z < zr; z++) {
        HRAD[z] = (((TTM[z + 1] - 2.0f * TTM[z] + TTM[z - 1])) * h2z * dzzdz[z] * dzzdz[z] + (TTM[z + 1] - TTM[z - 1]) * hz * d2zzdz2[z]) * RKAPA[z];
        HRAD[z] += (TTM[z + 1] - TTM[z - 1]) * hz * dzzdz[z] * DKAPA[z];
    }
    func_mul_extend<OPT::ADD, DIM::Z>(FORRANGE, FT, RO, HRAD, tmp); // do anonymous
#else
    tmp = OCV / REPR;
    func_d2_mul_array_z<OPT::ADD, DIM::X>(FORRANGE, FT, TT, RO, RKAPA, tmp); // do 519.1
    func_d2_mul_array_z<OPT::ADD, DIM::Y>(FORRANGE, FT, TT, RO, RKAPA, tmp); // do 519.2
    func_d2_mul_array_z<OPT::ADD, DIM::Z>(FORRANGE, FT, TT, RO, RKAPA, tmp); // do 519.3
    func_d1_mul_array_z<OPT::ADD, DIM::Z>(FORRANGE, FT, TT, RO, DKAPA, tmp); // do 519.4
#endif
    print_with_tag_index_all("after fluxes do 519", VIEW_INDX);
    if (RLAX > 0.0) {
        if ((TIMT > TSTART) && (TFAC > 1e-4f)) {
#if iTCon == 0 && izCon == 1
            if (myPE == nPE - 1) {
                if (TTM[zr - 1] <= TTM_MAX) {
                    TFAC -= DT / (3.0f * TOFF);
                } else {
                    TTM_MAX = TTM[zr - 1];
                    TFAC += DT / (3.0f * TOFF);
                    TFAC = min(TFAC, 2.0f);
                }
            }
            MPI_Bcast(&TFAC, 1, MPISIZE, nPE - 1, MPI_COMM_WORLD);
#elif iTCon == 1 && izCon == 0
            if (myPE == 0) {
                if (TTM[zl] >= TTM_MIN) {
                    TFAC -= DT / (3.0f * TOFF);
                } else {
                    TTM_MIN = TTM[zl];
                    TFAC += DT / (3.0f * TOFF);
                    TFAC = min(TFAC, 2.0f);
                }
            }
            MPI_Bcast(&TFAC, 1, MPISIZE, 0, MPI_COMM_WORLD);
#else
            exit_err("Relax: check boundary temperature");
#endif
            func_mul2<OPT::SET>(FORRANGE, WW2, RW, WW1, -2.5f);
            func_mul3<OPT::ADD>(FORRANGE, WW2, RW, UU, UU, -0.5f);
            func_mul3<OPT::ADD>(FORRANGE, WW2, RW, VV, VV, -0.5f);
            func_mul3<OPT::ADD>(FORRANGE, WW2, RW, WW, WW, -0.5f);
            horizontal_mean(FCONM, WW2);
            print_with_tag_index_all("Hmean ww2", VIEW_INDX);
            for (int z = zl; z < zr; z++) {
                FRADM[z] = (TTM[z + 1] - TTM[z - 1]) * hz * dzzdz[z] * RKAPA[z] / REPR;
            }
            if (myPEz == 0) {
                FRADM[zl] = (-3.0f * TTM[zl] + 4.0f * TTM[zl + 1] - TTM[zl + 2]) * hz * dzzdz[zl] * RKAPA[zl] / REPR;
            }
            if (myPEz == nPEz - 1) {
                FRADM[zr - 1] = (3.0f * TTM[zr - 1] - 4.0f * TTM[zr - 2] + TTM[zr - 3]) * hz * dzzdz[zr - 1] * RKAPA[zr - 1] / REPR;
            }
            int iTag = 100;
            if (myPEz == 0) {
                MPI_Status iStatus;
                MPI_Recv(&FRADM[nz - 1], 1, MPISIZE, myPE + nPEy, iTag, MPI_COMM_WORLD, &iStatus);
            } else if (myPEz == nPEz - 1) {
                MPI_Send(&FRADM[iz / 2], 1, MPISIZE, myPE - nPEy, iTag, MPI_COMM_WORLD);
            } else {
                MPI_Status iStatus;
                MPI_Sendrecv(&FRADM[iz / 2], 1, MPISIZE, myPE - nPEy, iTag, &FRADM[nz - 1], 1, MPISIZE, myPE + nPEy, iTag, MPI_COMM_WORLD, &iStatus);
            }
            double FTOT;
#if lRem
            FTOT = 8.07f * 0.4f * THETA;
#else
            FTOT = THETA / REPR;
#endif
            for (int z = iz / 2; z < nz; z++) ALPHA[z] = 1.0f - (FCONM[z] + FRADM[z]) / FTOT;
            for (int z = iz / 2; z < nz; z++) ALPHA[z] = ALPHA[z] / (1.0f + 4.0f * fabs(ALPHA[z]));
            CORRECT[iz / 2] = 0;
            for (int z = iz / 2 + 1; z < nz; z++) {
                CORRECT[z] = CORRECT[z - 1] + 0.5f * (ALPHA[z] + ALPHA[z - 1]) * (ZEE[z] - ZEE[z - 1]);
            }
            if (myPEy == 0) {
                ENDVAL[myPEz] = CORRECT[nz - 1];
                int iTag = 200;
                if (myPEz > 0) {
                    MPI_Send(&ENDVAL[myPEz], 1, MPISIZE, 0, iTag, MPI_COMM_WORLD);
                }
                if (myPEz == 0) {
                    ADDSUM[0] = ENDVAL[0];
                    for (int iPE = 1; iPE < nPEz; iPE++) {
                        MPI_Status iStatus;
                        MPI_Recv(&ADDSUM[iPE], 1, MPISIZE, iPE * nPEy, iTag, MPI_COMM_WORLD, &iStatus);
                    }
                }
            }
            MPI_Bcast(ADDSUM, nPEz, MPISIZE, 0, MPI_COMM_WORLD);
            print_with_tag_index_all("af fluex mpi_bcast", VIEW_INDX);
            double SUM = 0.0f;
// TODO : 检查下面代码是否有逻辑错误
#if iTCon == 0 && izCon == 1
            if (myPEz > 0) {
                for (int iPE = 0; iPE < myPEz; iPE++) {
                    SUM += ADDSUM[iPE];
                }
            }
#elif iTCon == 1 && izCon == 0
            for (int iPE = nPEz - 1; iPE >= myPEz; iPE--) {
                SUM -= ADDSUM[iPE];
                SUM -= ADDSUM[iPE];
            }
#else
            exit_err("Relax: check boundary temperature.");
#endif
            for (int z = zl; z < zr; z++) {
                CORRECT[z] += SUM;
            }
        } else {
            for (int z = zl; z < zr; z++) {
                CORRECT[z] = 0;
            }
        }
        func_base_extend<OPT::ADD, DIM::Z>(FORRANGE, FT, CORRECT, RLAX * TFAC);
        if (fmod(TIMT, 25.0) <= DT && (ICALL % 3 == 0)) {
#if iTCon == 0 && izCon == 1
            if (myPE == nPE - 1) {
                exit_err("yk not implement : sub 350");
                // write(*,'(A7,5(D15.6))') 'Relax: ',TIMT,TFAC,CORRECT(nz-iz/2),TTM(nz-iz/2),TTM_MAX
            }
#endif
#if iTCon == 1 && izCon == 0
            if (myPE == 0) {
                exit_err("yk not implement : sub 350");
                // write(*,'(A7,5(D15.6))') 'Relax: ',TIMT,TFAC,CORRECT(iz/2+1),TTM(iz/2+1),TTM_MIN
            }
#endif
        }
        ICALL++;
    }
#if iTCon == 2
    double CLN = -4.0f * log(2.0f) / HH / HH;
    int z, y, x;
    if (TC != 0.0) {
        FOR3D {
            FT[z][y][x] -= TP * OCV * 0.5f * (1.0f + tanh(log(2.0f) / QFH * (TIMT - TC))) * RO[z][y][x] *
                           exp(CLN * square(EXX[x] - XP) / HH) *
                           exp(CLN * square(WYY[y] - YP) / HH) *
                           exp(CLN * square(ZEE[z] - ZP) / HH);
        }
    } else {
        FOR3D {
            FT[z][y][x] -= TP * OCV * RO[z][y][x] *
                           exp(CLN * square(EXX[x] - XP) / HH) *
                           exp(CLN * square(WYY[y] - YP) / HH) *
                           exp(CLN * square(ZEE[z] - ZP) / HH);
        }
    }
#endif
    func_d2<OPT::ADD, DIM::X>(FORRANGE, FU, UU, ORE * C43); // do 530.1
    func_d2<OPT::ADD, DIM::Y>(FORRANGE, FU, UU, ORE);       // do 530.2
    func_d2<OPT::ADD, DIM::Z>(FORRANGE, FU, UU, ORE);       // do 530.3

    func_d2<OPT::ADD, DIM::X>(FORRANGE, FV, VV, ORE);       // do 540.1
    func_d2<OPT::ADD, DIM::Y>(FORRANGE, FV, VV, ORE * C43); // do 540.2
    func_d2<OPT::ADD, DIM::Z>(FORRANGE, FV, VV, ORE);       // do 540.3

    func_d2<OPT::ADD, DIM::X>(FORRANGE, FW, WW, ORE);       // do 550.1
    func_d2<OPT::ADD, DIM::Y>(FORRANGE, FW, WW, ORE);       // do 550.2
    func_d2<OPT::ADD, DIM::Z>(FORRANGE, FW, WW, ORE * C43); // do 550.3
    print_with_tag_index_all("af fluex do 550", VIEW_INDX);

    func_d1<OPT::SET, DIM::X>(FORRANGE_EXT_Y, WW1, UU); // do 650
    func_d1<OPT::SET, DIM::X>(FORRANGE_EXT_Y, WW2, VV); // do 670
#if !lShr
    func_d1<OPT::SET, DIM::X>(FORRANGE, WW3, WW); // do 690
    tmp = OCV * ORE;
    func_mul3<OPT::ADD>(FORRANGE, FT, RO, WW1, WW1, tmp * C43); // do 700.1
    func_mul3<OPT::ADD>(FORRANGE, FT, RO, WW2, WW2, tmp);       // do 700.2
    func_mul3<OPT::ADD>(FORRANGE, FT, RO, WW3, WW3, tmp);       // do 700.3
    func_mul2<OPT::SUB>(FORRANGE, FT, WW1, TT, OCV);            // do 710
#endif

    func_d1<OPT::ADD, DIM::Y>(FORRANGE, FU, WW2, C13 * ORE); // do 720
    func_d1<OPT::ADD, DIM::Y>(FORRANGE, FV, WW1, C13 * ORE); // do 730
    print_with_tag_index_all("af fluex do 730", VIEW_INDX);

#if !lShr
    func_d1<OPT::SET, DIM::Y>(FORRANGE, WW1, UU); // do 740
    tmp = OCV * ORE;
    func_mul3<OPT::ADD>(FORRANGE, FT, RO, WW1, WW1, tmp);     // do 750.1
    func_mul3<OPT::ADD>(FORRANGE, FT, RO, WW1, WW2, tmp * 2); // do 750.2
    print_with_tag_index_all("af before do 750", VIEW_INDX);
    print_with_tag_index_all("af fluex do 750", VIEW_INDX);
    func_d1<OPT::SET, DIM::Y>(FORRANGE, WW2, VV);             // do 755
    print_with_tag_index_all("af fluex do 755", VIEW_INDX);
    tmp = C43 * OCV * ORE;
    func_mul3<OPT::ADD>(FORRANGE, FT, RO, WW2, WW2, tmp); // do 760
    func_mul2<OPT::SUB>(FORRANGE, FT, WW2, TT, OCV);      // do 770
    func_d1<OPT::SET, DIM::X>(FORRANGE, WW1, UU);         // do 780
#endif

    func_d1<OPT::SET, DIM::Z>(FORRANGE_EXT_XY, WW3, WW); // do 790
    print_with_tag_index_all("af fluex do 790", VIEW_INDX);

#if !lShr
    tmp = C43 * OCV * ORE;
    func_mul3<OPT::ADD>(FORRANGE, FT, RO, WW3, WW3, tmp);  // do 820.1
    func_mul3<OPT::ADD>(FORRANGE, FT, RO, WW1, WW3, -tmp); // do 820.2
    func_mul3<OPT::ADD>(FORRANGE, FT, RO, WW1, WW2, -tmp); // do 820.3
    func_mul3<OPT::ADD>(FORRANGE, FT, RO, WW2, WW3, -tmp); // do 820.4
    func_mul2<OPT::SUB>(FORRANGE, FT, WW3, TT, OCV);       // do 830
#endif

    func_d1<OPT::ADD, DIM::X>(FORRANGE, FU, WW3, C13 * ORE); // do 840
    func_d1<OPT::ADD, DIM::Y>(FORRANGE, FV, WW3, C13 * ORE); // do 850

    func_d1<OPT::SET, DIM::Z>(FORRANGE_EXT_X, WW1, UU); // do 860
    func_d1<OPT::SET, DIM::Z>(FORRANGE_EXT_Y, WW2, VV); // do 880

#if !lShr
    func_d1<OPT::SET, DIM::Y>(FORRANGE, WW3, WW); // do 900
    tmp = OCV * ORE;
    func_mul3<OPT::ADD>(FORRANGE, FT, RO, WW1, WW1, tmp);        // do 910.1
    func_mul3<OPT::ADD>(FORRANGE, FT, RO, WW2, WW2, tmp);        // do 910.2
    func_mul3<OPT::ADD>(FORRANGE, FT, RO, WW3, WW3, tmp);        // do 910.3
    func_mul3<OPT::ADD>(FORRANGE, FT, RO, WW2, WW3, 2.0f * tmp); // do 910.4 // TODO: merge 910.2,3,4
#endif
    func_d1<OPT::ADD, DIM::X>(FORRANGE, FW, WW1, C13 * ORE); // do 920.1
    func_d1<OPT::ADD, DIM::Y>(FORRANGE, FW, WW2, C13 * ORE); // do 920.2
#if !lShr
    func_d1<OPT::SET, DIM::X>(FORRANGE, WW3, WW); // do 930
    tmp = OCV * ORE;
    func_mul3<OPT::ADD>(FORRANGE, FT, RO, WW1, WW3, 2.0f * tmp); // do 940
    print_with_tag_index_all("af fluex do 940", VIEW_INDX);
#endif

#if lRot
    func_base<OPT::ADD>(FORRANGE, FU, RV, OMZ);
    func_base<OPT::SUB>(FORRANGE, FV, RU, OMZ);
    func_base<OPT::ADD>(FORRANGE, FV, RW, OMX);
    func_base<OPT::SUB>(FORRANGE, FW, RV, OMX);
#endif

#if lMag
    func_d1<OPT::SET, DIM::X>(FORRANGE, RU, UU);
    func_d1<OPT::SET, DIM::X>(FORRANGE, RV, VV);
    func_d1<OPT::SET, DIM::X>(FORRANGE, RW, WW);

    func_mul2<OPT::SET>(FORRANGE, WW2, BX, RV);
    func_mul2<OPT::SUB>(FORRANGE, WW2, BY, RU);

    func_mul2<OPT::SET>(FORRANGE, WW3, BX, RW);
    func_mul2<OPT::SUB>(FORRANGE, WW3, BZ, RU);

    func_d1<OPT::SET, DIM::Y>(FORRANGE, RU, UU);
    func_d1<OPT::SET, DIM::Y>(FORRANGE, RV, VV);
    func_d1<OPT::SET, DIM::Y>(FORRANGE, RW, WW);

    func_mul2<OPT::SET>(FORRANGE, WW1, BY, RU);
    func_mul2<OPT::SUB>(FORRANGE, WW1, BX, RV);

    func_mul2<OPT::ADD>(FORRANGE, WW3, BY, RW);
    func_mul2<OPT::SUB>(FORRANGE, WW3, BZ, RV);

    func_d1<OPT::SET, DIM::Z>(FORRANGE, RU, UU);
    func_d1<OPT::SET, DIM::Z>(FORRANGE, RV, VV);
    func_d1<OPT::SET, DIM::Z>(FORRANGE, RW, WW);

    if (myPEz == 0) {
        for (int y = yl; y < yr; y++) {
            for (int x = xl; x < xr; x++) {
                RW[zl][y][x] = (4.0f * WW[zl + 1][y][x] - WW[zl + 2][y][x]) * hz * dzzdz[zl];
            }
        }
    } else if (myPEz == nPEz - 1) {
        for (int y = yl; y < yr; y++) {
            for (int x = xl; x < xr; x++) {
                RW[zr - 1][y][x] = (-4.0f * WW[zr - 2][y][x] + WW[zr - 3][y][x]) * hz * dzzdz[zr - 1];
            }
        }
    }

    func_mul2<OPT::ADD>(FORRANGE, WW1, BZ, RU);
    func_mul2<OPT::SUB>(FORRANGE, WW1, BX, RW);

    func_mul2<OPT::ADD>(FORRANGE, WW2, BZ, RV);
    func_mul2<OPT::SUB>(FORRANGE, WW2, BY, RW);

    func_d1<OPT::SET, DIM::X>(FORRANGE, RU, BX);
    func_d1<OPT::SET, DIM::X>(FORRANGE, RV, BY);
    func_d1<OPT::SET, DIM::X>(FORRANGE, RW, BZ);

    func_d1<OPT::SET, DIM::Y>(FORRANGE, TT, BX);
    func_mul2<OPT::ADD>(FORRANGE, FU, BY, TT, OBETA);
    func_mul2<OPT::SUB>(FORRANGE, FU, BY, RV, OBETA);
    func_mul2<OPT::SUB>(FORRANGE, FU, BZ, RW, OBETA);
    func_mul2<OPT::ADD>(FORRANGE, FV, BX, RV, OBETA);
    func_mul2<OPT::SUB>(FORRANGE, FV, BX, TT, OBETA);
    func_mul2<OPT::ADD>(FORRANGE, FW, BX, RW, OBETA);

    func_mul3<OPT::ADD>(FORRANGE, FT, RO, TT, TT, OBETA * ORM * OCV);
    func_mul3<OPT::ADD>(FORRANGE, FT, RO, RV, RV, OBETA * ORM * OCV);
    func_mul3<OPT::ADD>(FORRANGE, FT, RO, RW, RW, OBETA * ORM * OCV);
    func_mul3<OPT::ADD>(FORRANGE, FT, RO, RV, TT, -2.0f * OBETA * ORM * OCV);

    func_mul2<OPT::SUB>(FORRANGE, WW1, UU, RU);
    func_mul2<OPT::SUB>(FORRANGE, WW1, VV, TT);
    func_mul2<OPT::SUB>(FORRANGE, WW2, UU, RV);
    func_mul2<OPT::SUB>(FORRANGE, WW3, UU, RW);

    func_d2<OPT::SET, DIM::X>(FORRANGE, RU, BX);
    func_d2<OPT::SET, DIM::X>(FORRANGE, RV, BY);
    func_d2<OPT::SET, DIM::Y>(FORRANGE, TT, BX);

    func_base<OPT::ADD>(FORRANGE, WW1, RU, ORM);
    func_base<OPT::ADD>(FORRANGE, WW1, TT, ORM);
    func_base<OPT::ADD>(FORRANGE, WW2, RV, ORM);

    func_d2<OPT::SET, DIM::X>(FORRANGE, TT, BZ);

    func_base<OPT::ADD>(FORRANGE, WW3, TT, ORM);

    func_d1<OPT::SET, DIM::Z>(FORRANGE, RU, BX);
    func_d1<OPT::SET, DIM::Y>(FORRANGE, RV, BY);

    func_mul2<OPT::ADD>(FORRANGE, FU, BZ, RU, OBETA);
    func_mul2<OPT::SUB>(FORRANGE, FW, BX, RU, OBETA);

    func_mul3<OPT::ADD>(FORRANGE, FT, RO, RU, RU, OBETA * ORM * OCV);
    func_mul3<OPT::ADD>(FORRANGE, FT, RO, RU, RW, -2.0f * OBETA * ORM * OCV);

    func_mul2<OPT::SUB>(FORRANGE, WW1, WW, RU);
    func_mul2<OPT::SUB>(FORRANGE, WW2, VV, RV);

    func_d2<OPT::SET, DIM::Z>(FORRANGE, RU, BX);

    if (myPEz == 0) {
        for (int y = yl; y < yr; y++) {
            for (int x = xl; x < xr; x++) {
                RU[zl][y][x] = (-3.0f * BX[zl][y][x] + 4.0f * BX[zl + 1][y][x] - BX[zl + 2][y][x]) * hz * d2zzdz2[zl] +
                               (2.0f * BX[zl][y][x] - 5.0f * BX[zl + 1][y][x] + 4.0f * BX[zl + 2][y][x] - BX[zl + 3][y][x]) * h2z * dzzdz[zl] * dzzdz[zl];
            }
        }
    } else if (myPEz == nPEz - 1) {
        for (int y = yl; y < yr; y++) {
            for (int x = xl; x < xr; x++) {
                RU[zr - 1][y][x] = (3.0f * BX[zr - 1][y][x] - 4.0f * BX[zr - 2][y][x] + BX[zr - 3][y][x]) * hz * d2zzdz2[zr - 1] +
                                   (2.0f * BX[zr - 1][y][x] - 5.0f * BX[zr - 2][y][x] + 4.0f * BX[zr - 3][y][x] - BX[zr - 4][y][x]) * h2z * dzzdz[zr - 1] * dzzdz[zr - 1];
            }
        }
    }

    func_d2<OPT::SET, DIM::Y>(FORRANGE, RV, BY);
    func_base<OPT::ADD>(FORRANGE, WW1, RU, ORM);
    func_base<OPT::ADD>(FORRANGE, WW2, RV, ORM);
    func_d1<OPT::SET, DIM::Z>(FORRANGE, RV, BY);
    func_d1<OPT::SET, DIM::Y>(FORRANGE, RW, BZ);
    func_base<OPT::SET>(FORRANGE, RU, RV);
    func_base<OPT::SUB>(FORRANGE, RU, RW);
    func_mul2<OPT::ADD>(FORRANGE, FV, BZ, RU, OBETA);
    func_mul2<OPT::SUB>(FORRANGE, FW, BY, RU, OBETA);

    func_mul3<OPT::ADD>(FORRANGE, FT, RO, RV, RV, OBETA * ORM * OCV);
    func_mul3<OPT::ADD>(FORRANGE, FT, RO, RW, RW, OBETA * ORM * OCV);
    func_mul3<OPT::ADD>(FORRANGE, FT, RO, RV, RW, -2.0f * OBETA * ORM * OCV);

    func_mul2<OPT::SUB>(FORRANGE, WW2, WW, RV);
    func_mul2<OPT::SUB>(FORRANGE, WW3, VV, RW);
    func_d2<OPT::SET, DIM::Z>(FORRANGE, RV, BY);
    if (myPEz == 0) {
        for (int y = yl; y < yr; y++) {
            for (int x = xl; x < xr; x++) {
                RV[zl][y][x] = (-3.0f * BY[zl][y][x] + 4.0f * BY[zl + 1][y][x] - BY[zl + 2][y][x]) * hz * d2zzdz2[zl] +
                               (2.0f * BY[zl][y][x] - 5.0f * BY[zl + 1][y][x] + 4.0f * BY[zl + 2][y][x] - BY[zl + 3][y][x]) * h2z * dzzdz[zl] * dzzdz[zl];
            }
        }
    } else if (myPEz == nPEz - 1) {
        for (int y = yl; y < yr; y++) {
            for (int x = xl; x < xr; x++) {
                RV[zr - 1][y][x] = (3.0f * BY[zr - 1][y][x] - 4.0f * BY[zr - 2][y][x] + BY[zr - 3][y][x]) * hz * d2zzdz2[zr - 1] +
                                   (2.0f * BY[zr - 1][y][x] - 5.0f * BY[zr - 2][y][x] + 4.0f * BY[zr - 3][y][x] - BY[zr - 4][y][x]) * h2z * dzzdz[zr - 1] * dzzdz[zr - 1];
            }
        }
    }
    func_d2<OPT::SET, DIM::Y>(FORRANGE, RW, BZ);
    func_base<OPT::ADD>(FORRANGE, WW2, RV, ORM);
    func_base<OPT::ADD>(FORRANGE, WW3, RW, ORM);
    func_d1<OPT::SET, DIM::Z>(FORRANGE, RW, BZ);
    func_mul2<OPT::SUB>(FORRANGE, WW3, WW, RW);
    func_d2<OPT::SET, DIM::Z>(FORRANGE, RW, BZ);
    func_base<OPT::ADD>(FORRANGE, WW3, RW, ORM);
#endif
}

void horizontal_mean(double varM[nz], data var) {
    static double wwz[nz];
    int z, y, x;

#if nGrid == 0
    For(z, 0, nz) {
        wwz[z] = 0;
        for (y = iy / 2; y < ny - iy / 2; y++)
            for (x = ix / 2; x < nx - ix / 2; x++)
                wwz[z] += var[z][y][x];
        wwz[z] /= float(nPx * nPy);
    }
#else
    exit_err("Update spline interpolation");
#endif
    if (nPEy == 1) {
        for (int z = 0; z < nz; z++) varM[z] = wwz[z];
    } else {
        MPI_Allreduce(wwz, varM, nz, MPISIZE, MPI_SUM, MPI_COMM_WORLD);
    }
}

void BCON(FORRANGEPARAM) {
#if lShr
    double RRT = 0, RRB = 0, RPP = 0;
    RRT = 1.0f;
    RRB = -1.0f;
    RPP = 0.0f;
#endif
    if (myPEz == 0) {
#if lShr
        for (int y = yl; y < yr; y++) {
            for (int x = xl; x < xr; x++) {
                RU[zl][y][x] = 0.0f;
                RV[zl][y][x] = RRT * cos(RPP * TIMT) * RO[zl][y][x];
            }
        }
#else
        for (int y = yl; y < yr; y++) {
            for (int x = xl; x < xr; x++) {
                RU[zl][y][x] = RO[zl][y][x] * (C43 * RU[zl + 1][y][x] / RO[zl + 1][y][x] -
                                               C13 * RU[zl + 2][y][x] / RO[zl + 2][y][x]);
                RV[zl][y][x] = RO[zl][y][x] * (C43 * RV[zl + 1][y][x] / RO[zl + 1][y][x] -
                                               C13 * RV[zl + 2][y][x] / RO[zl + 2][y][x]);
            }
        }
#endif
        for (int y = yl; y < yr; y++) {
            for (int x = xl; x < xr; x++) {
                RW[zl][y][x] = 0.0f;
            }
        }
        if ((TP != 0.0) && (iTCon == 0 || iTCon == 1)) {
#if iTCon == 0 || iTCon == 1
            double TPR;
            if (TC != 0.0) {
                TPR = TP + (1.0f - TP) * exp(-TIMT / TC);
            } else {
                TPR = TP;
            }
    #if iTCon == 0
            double CLN = -4.0f * log(2.0f) / HH / HH;
            for (int y = yl; y < yr; y++) {
                for (int x = xl; x < xr; x++) {
                    TT[zl][y][x] = 1.0f - (1.0f - TPR) *
                                              exp(CLN * square(EXX[x] - XP)) / HH *
                                              exp(CLN * square(WYY[y] - YP)) / HH;
                }
            }
    #else
            double CLN = -4.0f * log(2.0f) / HH / HH;
            double DZ = 1.0f / double(nPz - 1) / dzzdz[zl];
            for (int y = yl; y < yr; y++) {
                for (int x = xl; x < xr; x++) {
                    TT[zl][y][x] = C43 * TT[zl + 1][y][x] -
                                   C13 * TT[zl + 2][y][x] -
                                   C23 * DZ * (THETA - (THETA - TPR) * exp(CLN * square(EXX[x] - XP)) / HH * exp(CLN * square(WYY[y] - YP)) / HH);
                }
            }
    #endif
#endif
        } else {
#if iTCon == 0 || iTCon == 2 || iTCon == 4
            for (int y = yl; y < yr; y++) {
                for (int x = xl; x < xr; x++) {
    #if lRem
                    TT[zl][y][x] = TU;
    #else
                    TT[zl][y][x] = 1.0f;
        #if iTCon == 4
                    if (x < nx / 2) {
                        TT[zl][y][x] = 1.0f - TP;
                    }
        #endif
    #endif
                }
            }
#endif
#if iTCon == 1
            double DZ = 1.0f / float(nPz - 1) / dzzdz[zl];
            for (int y = yl; y < yr; y++) {
                for (int x = xl; x < xr; x++) {
    #if lRem
                    TT[zl][y][x] = C43 * TT[zl + 1][y][x] -
                                   C13 * TT[zl + 2][y][x] + DZTU;
    #else
                    TT[zl][y][x] = C43 * TT[zl + 1][y][x] -
                                   C13 * TT[zl + 2][y][x] - C23 * DZ * THETA;
    #endif
                }
            }
#endif
        }
#if lMag && iBCon == 0
        for (int y = yl; y < yr; y++) {
            for (int x = xl; x < xr; x++) {
                BZ[zl][y][x] = 0.0f;
            }
        }
#endif
    }
    if (myPEz == nPEz - 1) {
#if lShr
        for (int y = yl; y < yr; y++) {
            for (int x = xl; x < xr; x++) {
                RU[zr - 1][y][x] = 0.0f;
                RV[zr - 1][y][x] = RRB * cos(RPP * TIMT) * RO[zr - 1][y][x];
            }
        }
#else
        for (int y = yl; y < yr; y++) {
            for (int x = xl; x < xr; x++) {
                RU[zr - 1][y][x] = RO[zr - 1][y][x] * (C43 * RU[zr - 2][y][x] / RO[zr - 2][y][x] -
                                                       C13 * RU[zr - 3][y][x] / RO[zr - 3][y][x]);
                RV[zr - 1][y][x] = RO[zr - 1][y][x] * (C43 * RV[zr - 2][y][x] / RO[zr - 2][y][x] -
                                                       C13 * RV[zr - 3][y][x] / RO[zr - 3][y][x]);
            }
        }
#endif
        for (int y = yl; y < yr; y++) {
            for (int x = xl; x < xr; x++) {
                RW[zr - 1][y][x] = 0.0f;
            }
        }
#if izCon == 0
        for (int y = yl; y < yr; y++) {
            for (int x = xl; x < xr; x++) {
                TT[zr - 1][y][x] = TB;
            }
        }
#elif izCon == 1
        double DZ = 1.0f / float(nPz - 1) / dzzdz[zr - 1];
        for (int y = yl; y < yr; y++) {
            for (int x = xl; x < xr; x++) {
    #if lRem
                TT[zr - 1][y][x] = C43 * TT[zr - 2][y][x] -
                                   C13 * TT[zr - 3][y][x] + DZTB;
    #else
                TT[zr - 1][y][x] = C43 * TT[zr - 2][y][x] -
                                   C13 * TT[zr - 3][y][x] -
                                   C23 * DZ * THETA / RKAPA[zr - 1];
    #endif
            }
        }
#else
        exit_err("BCON:  Invalid izCon");
#endif
#if lMag && iBCon == 0
        for (int y = yl; y < yr; y++) {
            for (int x = xl; x < xr; x++) {
                BZ[zr - 1][y][x] = 0.0f;
            }
        }
#endif
    }
}

#endif
