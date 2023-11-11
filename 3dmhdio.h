#ifndef _3DMHD_IO_H_
#define _3DMHD_IO_H_

#include "3dmhddata.h"
#include "3dmhdutils.h"

#include <cmath>
#include <fstream>
#include <iostream>

// local
char blanks[81], strng[81];
int nCase, nCaseP, nStart, nDump0, NSW;
double par1[96];
double PR, RY, ANG, BETA, POLYS;
double wMin[4], wMout[4];

double DT1, DT2, DT3, DT4;

void dumpData(std::string dataFile, int idx) {
    int rec = lMag ? 8 * idx : 5 * idx;
    std::ofstream fout(dataFile, std::ios::binary | std::ios::ate);
    fout.seekp(rec * iWord * nCore);
    save(fout, RU, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
    save(fout, RV, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
    save(fout, RW, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
    save(fout, TT, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
    save(fout, RO, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
#if lMag
    save(fout, BX, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
    save(fout, BY, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
    save(fout, BZ, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
#endif
    fout.close();
}

void dumpParam(std::string paramFile, int idx) {
    std::ofstream fout(paramFile, std::ios::binary | std::ios::ate);
    fout.seekp(idx * iWord * 96, std::ios::beg);
    fout.write(reinterpret_cast<char*>(par1), sizeof(double) * 96);
    fout.close();
}

void dump(std::string dataFile, std::string paramFile) {
    dumpData(dataFile, nDump0);
    if (myPE == 0) dumpParam(paramFile, nDump0);
    nDump0 += 1;
}

void loadParam(std::string paramFile, int idx) {
    std::ifstream fin(paramFile, std::ios::binary);
    fin.seekg(idx * iWord * 96, std::ios::beg);
    fin.read(reinterpret_cast<char*>(par1), sizeof(double) * 96);
    fin.close();
}

void loadData(std::string dataFile, int idx) {
    int rec = lMag ? 8 * idx : 5 * idx;
    std::ifstream fin(dataFile, std::ios::binary);
    fin.seekg(rec * iWord * nCore);
    load(fin, RU, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
    load(fin, RV, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
    load(fin, RW, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
    load(fin, TT, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
    load(fin, RO, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
#if lMag
    load(fin, BX, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
    load(fin, BY, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
    load(fin, BZ, iz / 2, nz - iz / 2, 1, ny - iy + 1, 1, nx - ix + 1);
#endif
    fin.close();
}

#define s4000(buffer)          sprintf(buffer, "------------------------------------------------------------------------------@")
#define s4001(buffer)          sprintf(buffer, "                                                                              @")
#define s4002(buffer)          sprintf(buffer, "%5sTHREE-D MAGNETOHYDRODYNAMICS@", "")
#define s4003(buffer, i)       sprintf(buffer, "%5sCase #%3d@", "", i)
#define s4004(buffer)          sprintf(buffer, "%5sStatic start@", "")
#define s4005(buffer)          sprintf(buffer, "%5sContinuation@", "")
#define s4006(buffer, i)       sprintf(buffer, "%5sPrevious case #%3d@", "", i)
#define s4007(buffer, x, y, z) sprintf(buffer, "%5sArray size: nx=%4d ny=%4d nz=%4d@", "", x, y, z)
#define s4008(buffer)          sprintf(buffer, "%5sRotating domain.@", "")
#define s4009(buffer)          sprintf(buffer, "%5sNonrotating domain.@", "")
#define s4010(buffer)          sprintf(buffer, "%5sMagnetic fields included.@", "")
#define s4011(buffer)          sprintf(buffer, "%5sNo magnetic fields.@", "")
#define s4012(buffer)          sprintf(buffer, "%5sHorizontally periodic.@", "")
#define s4013(buffer)          sprintf(buffer, "%5sStress free and impenetrable lower boundary.@", "")
#define s4014(buffer)          sprintf(buffer, "%5sConstant temperature lower boundary.@", "")
#define s4015(buffer)          sprintf(buffer, "%5sConstant thermal flux lower boundary.@", "")
#define s4016(buffer)          sprintf(buffer, "%5sConstant temperature upper boundary.@", "")
#define s4017(buffer)          sprintf(buffer, "%5sConstant thermal flux upper boundary.@", "")
#define s4059(buffer)          sprintf(buffer, "%5sEmbedded heat loss.@", "")
#define s4060(buffer)          sprintf(buffer, "%5sNo temperature perturbation.@", "")
#define s4061(buffer)          sprintf(buffer, "%5sEmbedded thermal.@", "")
#define s4062(buffer)          sprintf(buffer, "%5sDifferential upper boundary temperature.@", "")
#define s4018(buffer)          sprintf(buffer, "%5sParameters for this run:@", "")
#define s4019(buffer, i)       sprintf(buffer, "%10sRe    = %11.4E@", "", i)
#define s4020(buffer, i)       sprintf(buffer, "%10sPr    = %11.4E@", "", i)
#define s4021(buffer, i)       sprintf(buffer, "%10sTheta = %11.4E@", "", i)
#define s4022(buffer, i)       sprintf(buffer, "%10sGrav  = %11.4E@", "", i)
#define s4023(buffer, i)       sprintf(buffer, "%10sRo    = %11.4E@", "", i)
#define s4024(buffer, i)       sprintf(buffer, "%10sAngle = %11.4E@", "", i)
#define s4025(buffer, i)       sprintf(buffer, "%10sRm    = %11.4E@", "", i)
#define s4026(buffer, i)       sprintf(buffer, "%10sBeta  = %11.4E@", "", i)
#define s4027(buffer, i)       sprintf(buffer, "%10sZpen  = %11.4E@", "", i)
#define s4028(buffer, i)       sprintf(buffer, "%10sSigma = %11.4E@", "", i)
#define s4029(buffer, i)       sprintf(buffer, "%10sms    = %11.4E@", "", i)
#define s4030(buffer, i)       sprintf(buffer, "%10sTp    = %11.4E@", "", i)
#define s4031(buffer, i)       sprintf(buffer, "%10sXp    = %11.4E@", "", i)
#define s4032(buffer, i)       sprintf(buffer, "%10sYp    = %11.4E@", "", i)
#define s4089(buffer, i)       sprintf(buffer, "%10sZp    = %11.4E@", "", i)
#define s4033(buffer, i)       sprintf(buffer, "%10sXmax  = %11.4E@", "", i)
#define s4034(buffer, i)       sprintf(buffer, "%10sYmax  = %11.4E@", "", i)
#define s4035(buffer, i)       sprintf(buffer, "%10sZmax  = %11.4E@", "", i)
#define s4036(buffer)          sprintf(buffer, "%5sTime stepping safety factor:@", "")
#define s4037(buffer, i)       sprintf(buffer, "%10sSFF   = %11.4E@", "", i)
#define s4038(buffer)          sprintf(buffer, "%5sCalculation started from a static state.@", "")
#define s4039(buffer, i)       sprintf(buffer, "%5sInitial temperature perturbations: %11.4E@", "", i)
#define s4040(buffer, i)       sprintf(buffer, "%5sInitial magnetic layer amplitude: %11.4E@", "", i)
#define s4041(buffer, i, s)    sprintf(buffer, "%5sStarted from dump %d of file %20s@", "", i, s)
#define s4042(buffer, i)       sprintf(buffer, "%5sPrevious simulation time: %10.4E@", "", i)
#define s4043(buffer, x, y, z) sprintf(buffer, "%5sPrevious size: nx=%4d ny=%4d nz=%4d@", "", x, y, z)
#define s4044(buffer, s)       sprintf(buffer, "%5sThe results are stored in file %15s@", "", s)
#define s4045(buffer, i)       sprintf(buffer, "%5sTotal number of iterations required = %8d@", "", i)
#define s4046(buffer, i)       sprintf(buffer, "%5sVariables dumped every %6d iterations.@", "", i)
#define s4047(buffer, i)       sprintf(buffer, "%5sAdvected fast-mode time: %10.4E@", "", i)
#define s4048(buffer, i)       sprintf(buffer, "%5sThermal diffusion time:  %10.4E@", "", i)
#define s4049(buffer, i)       sprintf(buffer, "%5sViscous diffusion time:  %10.4E@", "", i)
#define s4050(buffer, i)       sprintf(buffer, "%5sMagnetic diffusion time: %10.4E@", "", i)
#define s4051(buffer, i)       sprintf(buffer, "%5sInitial time-step value: %10.4E@", "", i)
#define s4052(buffer, i)       sprintf(buffer, "%5sIterations per crossing time: %8d@", "", i)
#define s4053(buffer, i)       sprintf(buffer, "%5sIteration %7d completed@", "", i)
#define s4054(buffer)          sprintf(buffer, "%5s-----------------------------------@", "")
#define s4055(buffer, i)       sprintf(buffer, "%5sTotal simulation time:   %10.4E@", "", i)
#define s4056(buffer, i)       sprintf(buffer, "%5sPresent simulation time: %10.4E@", "", i)
#define s4057(buffer, i)       sprintf(buffer, "%5sMaximum Mach number:     %10.4E@", "", i)
#define s4058(buffer, i)       sprintf(buffer, "%5sPresent time step:       %10.4E@", "", i)
#define s4090(buffer)          sprintf(buffer, "%15sSIMULATION COMPLETE@", "")
#define s4100(buffer, i, j)    sprintf(buffer, "%5sTime step %11.4E  Shut down at step %8d@", "", i, j)

void logGreet(std::string logFile) {
    std::ofstream outLog(logFile, std::ios::out);
    if (!outLog.is_open()) exit_err("Can't open .lis");
    s4000(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4001(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4001(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4002(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4003(strng, nCase), standard(strng), outLog.write(strng, sizeof(char) * 80);
    if (nStart == 0) {
        s4004(strng);
    } else if (nCase == nCaseP) {
        s4005(strng);
    } else {
        s4006(strng, nCaseP);
    }
    standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4001(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4001(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4000(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4001(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4007(strng, nPx, nPy, nPz), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4001(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    lRot ? s4008(strng) : s4009(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    lMag ? s4010(strng) : s4011(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4012(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4013(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    izCon ? s4014(strng) : s4015(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    switch (iTCon) {
        case 0: s4016(strng); break;
        case 1: s4017(strng); break;
        case 2: s4059(strng); break;
        case 3: s4061(strng); break;
        case 4: s4062(strng); break;
        default: s4060(strng); break;
    }
    standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4001(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4018(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4019(strng, RE), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4020(strng, PR), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4021(strng, THETA), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4022(strng, GRAV), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4023(strng, RY), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4024(strng, ANG), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4025(strng, RM), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4026(strng, BETA), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4027(strng, PZP), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4028(strng, SIGMA), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4029(strng, POLYS), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4030(strng, TP), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4031(strng, XP), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4032(strng, YP), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4089(strng, ZP), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4033(strng, xMax), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4034(strng, yMax), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4035(strng, zMax), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4001(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4036(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4037(strng, SFF), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4001(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    if (nStart == 0) {
        s4038(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
        s4039(strng, AMPT), standard(strng), outLog.write(strng, sizeof(char) * 80);
        lMag ? s4040(strng, AMPB) : s4011(strng);
        standard(strng), outLog.write(strng, sizeof(char) * 80);
    } else {
        s4041(strng, nStart, fInp.c_str()), standard(strng), outLog.write(strng, sizeof(char) * 80);
        s4042(strng, TIMI), standard(strng), outLog.write(strng, sizeof(char) * 80);
        s4043(strng, nPx, nRy * nPEy, nRz * nPEz), standard(strng), outLog.write(strng, sizeof(char) * 80);
    }
    s4001(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4044(strng, fOut.c_str()), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4045(strng, nTotal), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4046(strng, nStep0), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4001(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4047(strng, DT1), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4048(strng, DT2), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4049(strng, DT3), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4050(strng, DT4), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4051(strng, DT), standard(strng), outLog.write(strng, sizeof(char) * 80);
    int nSound = round(1.0f / DT);
    s4052(strng, nSound), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4001(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4000(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4000(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    outLog.close();
}

void logEpoch(std::string logFile) {
    std::ofstream outLog(logFile, std::ios::out | std::ios::app);
    s4001(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4053(strng, nIt), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4054(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4055(strng, TIMT), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4056(strng, TIMC), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4057(strng, UMACH), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4058(strng, DT), standard(strng), outLog.write(strng, sizeof(char) * 80);
    outLog.close();
}

void logFatal(std::string logFile) {
    std::ofstream outLog(logFile, std::ios::out | std::ios::app);
    s4000(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4100(strng, DT, nIt - 1), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4000(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    s4000(strng), standard(strng), outLog.write(strng, sizeof(char) * 80);
    outLog.close();
}

#endif // _3DMHD_IO_H_
