#ifndef _3DMHD_PHISIC_H_
#define _3DMHD_PHISIC_H_

#include "3dmhddata.h"
#include "3dmhdparam.h"
#include "3dmhdset.h"
#include "3dmhdutils.h"

#include <cstring>
#include <fstream>
#include <math.h>
#include <mpi.h>
#include <vector>

#define ICOORD_ERROR "mkGrid: Invalid ICOORD number"
#define NGRID_ERROR  "mkGrid: Invalid nGrid number"

void trapzd(double (*func)(double), double a, double b, double& s, int n);
void polint(double xa[], double ya[], int n, double x, double& y, double& dy);
void mmid(std::vector<double>& Y, const std::vector<double>& DYDX, int NVAR, double XS, double HTOT, int NSTP, double* YOUT);
void rzextr(int IEST, double XEST, double* YEST, double* YZ, std::vector<std::vector<double>>& d, double* DY, int NV, int NUSE);
void derivs(double ZZ, double* Y, int NV, double* DYDX);
void ODEINT(double* YSTART, double X1, double X2);
void RKQC(double* Y, double* DYDX, int N, double& X, double HTRY, double EPS, double* YSCAL, double& HDID, double& HNEXT);
void RK4(double* Y, double* DYDX, int N, double X, double H, double* YOUT);
void tube();

defineDimSwitch(const double&, uMax, xMax, yMax, zMax);
defineDimSwitch(const double&, uu1, XX1, YY1, ZZ1);
defineDimSwitch(const double&, uu2, XX2, YY2, ZZ2);
defineDimSwitch(const double&, ua, XA, YA, ZA);
defineDimSwitch(const double&, ub, XB, YB, ZB);
defineDimSwitch(const double&, uc, XC, YC, ZC);
defineDimSwitch(const double&, ud, XD, YD, ZD);

template <DIM dim> void mkGrid(double sCode, double& sPhys, double& dssds, double& d2ssds2) {
#if nGrid == 0
    double sMax = uMax<dim>();
    sPhys = sCode * sMax;
    dssds = 1.0f / sMax;
    d2ssds2 = 0.0f;
#elif nGrid == 1
    double sMax = uMax<dim>();
    double s1 = uu1<dim>(), s2 = uu2<dim>();
    double a1 = (s2 - s1) / sMax;
    double a = atan(s1);
    double a2 = atan(s2) - a3;
    if (sCode == 0.0) {
        sPhys = 0.0f;
    } else if (sCode == 1.0) {
        sPhys = sMax;
    } else {
        sPhys = (tan(a2 * sCode + a3) - s1) / a1;
    }
    double denominator = 1.0f + (s1 + a1 * sPhys) * (s1 + a1 * sPhys);
    dssds = a1 / a2 / denominator;
    d2ssds2 = -2.0f * a1 * a1 * (s1 + a1 * sPhys) / (a2 * denominator * denominator);
#elif nGrid == 2
    double a = ua<dim>(), b = ub<dim>(), c = uc<dim>(), d = ud<dim>();
    if (c <= 1.0E-09) {
        sPhys = sCode * sMax;
        dssds = 1.0f / sMax;
        d2ssds2 = 0.0f;
    } else {
        double sH = (1.0f - d) * PI / (d * atan((a - b) * c) + 2.0f * atan(b * c) - d * atan((a + b) * c));
        double sK = 2.0f * c * PI * sMax / (2.0f * c * PI + sH * (2.0f * c * (atan((a - b) * c) * (a - b) - atan((a + b) * c) * (a + b) + atan((a - b - 1.0f) * c) + atan((a + b - 1.0f) * c) * (a + b - 1.0f) + atan((1.0f - a + b) * c) * (a - b)) - log(1.0f + (a - b) * (a - b) * c * c) + log(1.0f + (a + b) * (a + b) * c * c) - log(1.0f + (a + b - 1.0f) * (a + b - 1.0f) * c * c) + log(1.0f + (1.0f - a + b) * (1.0f - a + b) * c * c)));
        dssds = 1.0f / (sK * (1.0f + sH / PI * (atan(c * (sCode - a - b)) + atan(c * (a - sCode - b)))));
        d2ssds2 = -dssds * dssds * dssds * c * sH * sK / PI * (1.0f / (1.0f + c * c * (a + b - sCode) * (a + b - sCode)) - 1.0f / (1.0f + c * c * (a - b - sCode) * (a - b - sCode)));
        sPhys = sK / (2.0f * c * PI) * (2.0f * c * PI * sCode + sH * (2.0f * c * (atan((a - b) * c) * (a - b) - atan((a + b) * c) * (a + b) + atan((a - b - sCode) * c) * sCode + atan((a + b - sCode) * c) * (a + b - sCode) + atan((sCode - a + b) * c) * (a - b)) - log(1.0f + (a - b) * (a - b) * c * c) + log(1.0f + (a + b) * (a + b) * c * c) - log(1.0f + (a + b - sCode) * (a + b - sCode) * c * c) + log(1.0f + (sCode - a + b) * (sCode - a + b) * c * c)));
    }
    if (sCode == 0.0) sPhys = 0.0f;
    if (sCode == 1.0f) sPhys = sMax;
    break;
#else
    exit_err(NGRID_ERROR);
#endif
}

void xJacobi() {
    // Returns the x coordinate Jacobian transformation elements.
    const double ORX = 1.0f / float(nx - ix - 1);
    for (int i = 0; i < nx - ix; ++i) EXX[i + 1] = float(i) * ORX;
    for (int i = 1; i < nx - ix + 1; ++i) {
        double sPhys, dssds, d2ssds2;
        mkGrid<DIM::X>(EXX[i], sPhys, dssds, d2ssds2);
        EXX[i] = sPhys;
        dxxdx[i] = dssds;
        d2xxdx2[i] = d2ssds2;
        DDX[i] = ORX * (1.0f / dssds);
    }
}

void yJacobi() {
    // Full computational grid with evenly spaced y values between zero and one
    const double ORY = 1.0f / float(nPy - 1);
    static double YY[nPy];
    for (int j = 0; j < nPy; ++j) YY[j] = float(j) * ORY;
    double sPhys, dssds, d2ssds2;
    if (myPEy == 0) {
        for (int j = 0; j < iy / 2; ++j) {
            mkGrid<DIM::Y>(0.0f, sPhys, dssds, d2ssds2);
            WYY[j] = sPhys, dyydy[j] = dssds, d2yydy2[j] = d2ssds2;
        }
        for (int j = iy / 2; j < nRy + iy; ++j) {
            mkGrid<DIM::Y>(YY[j - iy / 2], sPhys, dssds, d2ssds2);
            WYY[j] = sPhys, dyydy[j] = dssds, d2yydy2[j] = d2ssds2;
        }
    } else if (myPEy == nPEy - 1) {
        for (int j = 0; j < nRy + iy / 2; ++j) {
            mkGrid<DIM::Y>(YY[j + nPy - nRy - iy / 2], sPhys, dssds, d2ssds2);
            WYY[j] = sPhys, dyydy[j] = dssds, d2yydy2[j] = d2ssds2;
        }
        for (int j = nRy + iy / 2; j < nRy + iy; ++j) {
            mkGrid<DIM::Y>(1.0f, sPhys, dssds, d2ssds2);
            WYY[j] = sPhys, dyydy[j] = dssds, d2yydy2[j] = d2ssds2;
        }
    } else {
        for (int j = 0; j < ny; ++j) {
            mkGrid<DIM::Y>(YY[j + myPEy * nRy - iy / 2], sPhys, dssds, d2ssds2);
            WYY[j] = sPhys, dyydy[j] = dssds, d2yydy2[j] = d2ssds2;
        }
    }
    for (int j = 0; j < ny; ++j) DDY[j] = ORY * (1.0f / dyydy[j]);
}

void zJacobi() {
    const double ORZ = 1.0f / float(nPz - 1);
    static double ZZ[nPz];
    for (int k = 0; k < nPz; ++k) ZZ[k] = float(k) * ORZ;
    double sPhys, dssds, d2ssds2;
    if (myPEz == 0) {
        for (int k = 0; k < iz / 2; ++k) {
            mkGrid<DIM::Z>(0.0f, sPhys, dssds, d2ssds2);
            ZEE[k] = sPhys, dzzdz[k] = dssds, d2zzdz2[k] = d2ssds2;
        }
        for (int k = iz / 2; k < nRz + iz; ++k) {
            mkGrid<DIM::Z>(ZZ[k - iz / 2], sPhys, dssds, d2ssds2);
            ZEE[k] = sPhys, dzzdz[k] = dssds, d2zzdz2[k] = d2ssds2;
        }
    } else if (myPEz == nPEz - 1) {
        for (int k = 0; k < nRz + iz / 2; ++k) {
            mkGrid<DIM::Z>(ZZ[k + nPz - nRz - iz / 2], sPhys, dssds, d2ssds2);
            ZEE[k] = sPhys, dzzdz[k] = dssds, d2zzdz2[k] = d2ssds2;
        }
        for (int k = nRz + iz / 2; k < nRz + iz; ++k) {
            mkGrid<DIM::Z>(1.0f, sPhys, dssds, d2ssds2);
            ZEE[k] = sPhys, dzzdz[k] = dssds, d2zzdz2[k] = d2ssds2;
        }
    } else {
        for (int k = 0; k < nz; ++k) {
            mkGrid<DIM::Z>(ZZ[k + myPEz * nRz - iz / 2], sPhys, dssds, d2ssds2);
            ZEE[k] = sPhys, dzzdz[k] = dssds, d2zzdz2[k] = d2ssds2;
        }
    }
    for (int k = 0; k < nz; ++k) DDZ[k] = ORZ * (1.0f / dzzdz[k]);
}

double FF(double x) {
    return A * A * (x * x * x * x * x) * square(std::exp(-x * x) - C) /
           square(A * (x * x * x) + 1.0f);
}

void qromb(double (*func)(double), double a, double b, double& ss) {
    const int JMAX = 20;
    const int JMAXP = JMAX + 1;
    const int k = 5;
    const int KM = k - 1;
    const double EPS = 1.0e-9f;
    static double h[JMAXP], s[JMAXP];
    h[0] = 1.0f;
    for (int j = 0; j < JMAX; j++) {
        trapzd(func, a, b, s[j], j);
        if (j >= KM) {
            double dss;
            polint(&h[j - KM], &s[j - KM], k, 0.0f, ss, dss);
            if (fabs(dss) <= EPS * fabs(ss)) return;
        }
        s[j + 1] = s[j];
        h[j + 1] = 0.25f * h[j];
    }
    exit_err("Too many steps in qromb");
}

void polint(double xa[], double ya[], int n, double x, double& y, double& dy) {
    const int NMAX = 10;
    int ns = 0;
    double dif, dift;
    static double c[NMAX], d[NMAX];

    dif = std::abs(x - xa[0]);
    for (int i = 0; i < n; i++) {
        double dift = std::abs(x - xa[i]);
        if (dift < dif) ns = i, dif = dift;
        c[i] = d[i] = ya[i];
    }
    y = ya[ns--];
    for (int m = 0; m < n; m++) {
        for (int i = 0; i < n - m; i++) {
            double ho = xa[i] - x;
            double hp = xa[i + m] - x;
            double w = c[i + 1] - d[i];
            double den = ho - hp;
            if (den == 0.0) exit_err("Failure in polint");
            den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
        }
        y += (dy = (2 * ns < (n - m) ? c[ns + 1] : d[ns--]));
    }
}

void trapzd(double (*func)(double), double a, double b, double& s, int n) {
    if (n == 0) {
        s = 0.5f * (b - a) * (func(a) + func(b));
    } else {
        int it = std::pow(2, n - 1);
        double tnm = float(it);
        double del = (b - a) / tnm;
        double x = a + 0.5f * del;
        double sum = 0.0f;
        for (int j = 0; j < it; j++) {
            sum += func(x);
            x += del;
        }
        s = 0.5f * (s + (b - a) * sum / tnm);
    }
}

// CHECKED
void bsstep(double* Y, double* DYDX, int NV, double& X, double HTRY, double EPS, double& HDID, double& HNEXT) {
    const int NMAX = 10;
    const int IMAX = 11;
    const int NUSE = 7;
    const double ONE = 1.0f;
    const double SHRINK = 0.95f;
    const double GROW = 1.2f;

    double H = HTRY;
    double XSAV = X;

    static double YERR[NMAX], YSeq[NMAX];

    std::vector<double> YSAV(NMAX), DYSAV(NMAX);
    std::vector<int> NSeq{2, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96};
    std::vector<std::vector<double>> D(NMAX, std::vector<double>(NUSE));

    for (int i = 0; i < NV; ++i) YSAV[i] = Y[i], DYSAV[i] = DYDX[i];

    while (1) {
        for (int i = 0; i < IMAX; ++i) {
            mmid(YSAV, DYSAV, NV, XSAV, H, NSeq[i], YSeq);
            double XEST = square(H / NSeq[i]);
            rzextr(i, XEST, YSeq, Y, D, YERR, NV, NUSE);

            double ERRMAX = 0.0f;
            for (int j = 0; j < NV; ++j) {
                ERRMAX = std::max(ERRMAX, std::abs(YERR[j] / Y[j]));
            }
            ERRMAX /= EPS;
            if (ERRMAX < 1) {
                X = X + H;
                HDID = H;
                if (i == NUSE) {
                    HNEXT = H * SHRINK;
                } else if (i == NUSE - 1) {
                    HNEXT = H * GROW;
                } else {
                    HNEXT = (H * NSeq[NUSE - 1]) / NSeq[i];
                }
                return;
            }
        }
        H = 0.25f * H / std::pow(2.0f, (IMAX - NUSE) / 2);
        if (X + H == X) exit_err("BSSTEP: Step size underflow");
    }
}

// CHECKED
void mmid(std::vector<double>& Y, const std::vector<double>& DYDX, int NVAR, double XS, double HTOT, int NSTP, double* YOUT) {
    const int NMAX = 10;
    static double YM[NMAX], YN[NMAX];

    double H = HTOT / NSTP, X = XS + H;
    for (int i = 0; i < NVAR; ++i) YM[i] = Y[i], YN[i] = Y[i] + H * DYDX[i];

    // 假设 derivs 函数是外部定义的函数，根据你的实际情况调用它
    derivs(X, YN, NVAR, YOUT);
    double H2 = 2.0f * H;
    for (int N = 2; N <= NSTP; ++N) {
        for (int i = 0; i < NVAR; ++i) {
            double SWAP = YM[i] + H2 * YOUT[i];
            YM[i] = YN[i];
            YN[i] = SWAP;
        }
        X = X + H;
        derivs(X, YN, NVAR, YOUT);
    }
    for (int i = 0; i < NVAR; ++i) YOUT[i] = 0.5f * (YM[i] + YN[i] + H * YOUT[i]);
}

// CHECKED
void rzextr(int IEST, double XEST, double* YEST, double* YZ, std::vector<std::vector<double>>& D, double* DY, int NV, int NUSE) {
    const int IMAX = 11;
    const int NMAX = 10;
    const int NCOL = 7;

    // X 需要在多次调用中保存值  fortran : save x
    static double X[IMAX], FX[NCOL];

    X[IEST] = XEST;
    if (IEST == 0) {
        for (int j = 0; j < NV; ++j) YZ[j] = D[j][0] = DY[j] = YEST[j];
    } else {
        int M1 = std::min(IEST + 1, NUSE);
        for (int k = 1; k < M1; ++k) {
            FX[k] = X[IEST - k] / XEST;
        }

        for (int j = 0; j < NV; ++j) {
            double YY = YEST[j];
            double V = D[j][0];
            double C = YY;
            double DDY;
            D[j][0] = YY;

            for (int k = 1; k < M1; ++k) {
                double B1 = FX[k] * V;
                double B = B1 - C;
                if (B != 0.0) {
                    B = (C - V) / B;
                    DDY = C * B;
                    C = B1 * B;
                } else {
                    DDY = V;
                }
                V = D[j][k];
                D[j][k] = DDY;
                YY = YY + DDY;
            }
            DY[j] = DDY;
            YZ[j] = YY;
        }
    }
}

// CHECKED
void derivs(double zz, double* Y, int NV, double* DYDX) {
    double z, dzzdz, d2zzdz2;
    mkGrid<DIM::Z>(zz, z, dzzdz, d2zzdz2);
    double RKAP = 1.0f;
    if (PZP != 0.0) RKAP += (RKAPST - 1.0f) / 2.0f * (1.0f + tanh((zz - PZP) / SIGMA));
    // if(myPE==0) printf("RKAP=%.20lf\n",RKAP);
    DYDX[0] = THETA / RKAP / dzzdz;
    DYDX[1] = GRAV / (1.0f + Y[0]) / dzzdz;
}

void getbackground(double x, double& t, double& rho) {
    static double YSTART[2];
    YSTART[0] = YSTART[1] = 1;
    ODEINT(YSTART, 2.0f, x);
    t = YSTART[0], rho = YSTART[1] / YSTART[0];
}

/**
 * this function is defined in 3dmhdset.f, original version:
 *
c========================================================================
c  Specification of radiative conductivity (M. Rempel).
c========================================================================
      subroutine kappa(z,rho,t,krad)
      implicit none

      real*8 z,rho,t,alpha,beta,krad

c      beta=0.0
c      krad=t**6*sqrt(t)/rho**2
c      krad=(beta+krad)/(beta+1.0f)
c      krad=krad+(1.25-krad)*exp(-(z/0.25)**2)

      krad=rho
c      krad=(0.9**10+krad**10)**0.1

      krad=krad+(1.25-krad)*exp(-5.0*z)

      return
      end
c========================================================================
 *
*/
double kappa(double z, double rho, double t) {
    return rho + (1.25f - rho) * exp(-5.0f * z);
}

void derivs1(double x, double yt[2], double dyt[2]) {
    const double nad = 0.41f;
    const double os = 1.0f;
    double krad = kappa(x, yt[1] / yt[0], yt[0]);
    double nrad = 0.4f / krad;
    if (nrad <= os * nad) {
        dyt[0] = nrad;
        dyt[1] = yt[1] / yt[0];
    } else {
        dyt[0] = nad;
        dyt[1] = yt[1] / yt[0];
    }
}

// void (*derivs1)(double, double[], double[])
// void (*RKQC)(double[], double[], int, double&, double&, double, double[], double&, double&)
void ODEINT(double* YSTART, double X1, double X2) {
    const int NVAR = 2, NMAX = 2, MAXSTP = 1000;
    const int KMAX = 1;
    const double TINY = 1.0e-30f;
    const double DXSAV = 1e-4; // literal as double
    const double EPS = 1e-8;   // literal as double
    const double H1 = 1e-4;    // literal as double
    const double HMIN = 1e-50; // literal as double
    int nok = 0, nbad = 0, KOUNT = 0;

    static double YSCAL[NMAX], Y[NMAX], DYDX[NMAX], XP[100], YP[NVAR][100];
    double X = X1, H = (X2 > X1) ? H1 : -H1;

    for (int i = 0; i < NVAR; ++i) Y[i] = YSTART[i];
    double XSAV = X - DXSAV * 2.0f;
    for (int NSTP = 1; NSTP <= MAXSTP; ++NSTP) {
        derivs1(X, Y, DYDX);
        for (int i = 0; i < NVAR; ++i) {
            YSCAL[i] = std::abs(Y[i]) + std::abs(H * DYDX[i]) + TINY;
        }
        if (KOUNT < KMAX - 1 && std::abs(X - XSAV) > std::abs(DXSAV)) {
            XSAV = X;
            XP[KOUNT] = X;
            for (int i = 0; i < NVAR; ++i) YP[i][KOUNT] = Y[i];
            KOUNT++;
        }
        if ((X + H - X2) * (X + H - X1) > 0) H = X2 - X;
        double HDID, HNEXT;
        RKQC(Y, DYDX, NVAR, X, H, EPS, YSCAL, HDID, HNEXT);
        if (HDID == H) {
            nok++;
        } else {
            nbad++;
        }
        if ((X - X2) * (X2 - X1) >= 0) {
            // for (int i = 0; i < NVAR; ++i) YSTART[i] = Y[i];
            if (KMAX != 0) {
                XP[KOUNT] = X;
                for (int i = 0; i < NVAR; ++i) YP[i][KOUNT] = Y[i];
                KOUNT++;
            }
            // 将需要用到的数据返回 getbackground
            for (int i = 0; i < NVAR; ++i) YSTART[i] = YP[i][KOUNT - 1];
            return;
        }
        if (std::abs(HNEXT) < HMIN) exit_err("Stepsize smaller than minimum.");
        H = HNEXT;
    }
    exit_err("Too many steps.");
}

void RKQC(double* Y, double* DYDX, int N, double& X, double HTRY, double EPS, double* YSCAL, double& HDID, double& HNEXT) {
    const int MAXN = 10;
    const double FCOR = 1.0f / 15.0f;
    const double ONE = 1.0f;
    const double SAFETY = 0.9f;
    const double ERRCON = 6.0e-4f;
    static double YTEMP[MAXN], YSAV[MAXN], DYSAV[MAXN];

    const double PGROW = -0.20f;
    const double PSHRNK = -0.25f;
    double XSAV = X;

    for (int i = 0; i < N; i++) YSAV[i] = Y[i], DYSAV[i] = DYDX[i];

    double H = HTRY;

    while (true) {
        double HH = 0.5f * H;
        RK4(YSAV, DYSAV, N, XSAV, HH, YTEMP);
        X = XSAV + HH;
        derivs1(X, YTEMP, DYDX);
        RK4(YTEMP, DYDX, N, X, HH, Y);
        X = XSAV + H;
        if (X == XSAV) {
            exit_err("Stepsize not significant in RKQC.");
            break;
        }
        RK4(YSAV, DYSAV, N, XSAV, H, YTEMP);
        double ERRMAX = 0.0f;
        for (int i = 0; i < N; i++) {
            YTEMP[i] = Y[i] - YTEMP[i];
            ERRMAX = std::max(ERRMAX, fabs(YTEMP[i] / YSCAL[i]));
        }
        ERRMAX = ERRMAX / EPS;
        if (ERRMAX > 1.0) {
            H = SAFETY * H * pow(ERRMAX, PSHRNK);
        } else {
            HDID = H;
            if (ERRMAX > ERRCON) {
                HNEXT = SAFETY * H * pow(ERRMAX, PGROW);
            } else {
                HNEXT = 4.0f * H;
            }
            break;
        }
    }
    for (int i = 0; i < N; i++) Y[i] += YTEMP[i] * FCOR;
}

void RK4(double* Y, double* DYDX, int N, double X, double H, double* YOUT) {
    const int NMAX = 10;
    static double YT[NMAX], DYT[NMAX], DYM[NMAX];
    double HH = H * 0.5f, H6 = H / 6.0f, XH = X + HH;
    for (int i = 0; i < N; i++) YT[i] = Y[i] + HH * DYDX[i];
    derivs1(XH, YT, DYT);
    for (int i = 0; i < N; i++) YT[i] = Y[i] + HH * DYT[i];
    derivs1(XH, YT, DYM);
    for (int i = 0; i < N; i++) YT[i] = Y[i] + H * DYM[i], DYM[i] = DYT[i] + DYM[i];
    derivs1(X + H, YT, DYT);
    for (int i = 0; i < N; i++) YOUT[i] = Y[i] + H6 * (DYDX[i] + DYT[i] + 2.0f * DYM[i]);
}

double RAN2(int& IDUM) {
    const int M = 714025, IA = 1366, IC = 150889;
    const double RM = 1.0f / M;
    static int IIY = 0;
    static int IIR[97];

    if (IDUM < 0) {
        IDUM = fmod(IC - IDUM, M);
        for (int j = 0; j < 97; ++j) {
            IDUM = fmod(IA * IDUM + IC, M);
            IIR[j] = IDUM;
        }
        IDUM = fmod(IA * IDUM + IC, M);
        IIY = IDUM;
    }

    int j = (97 * IIY) / M;
    if (j >= 97 || j < 0) exit_err("RAN2: Error.");

    IIY = IIR[j];
    double randomValue = float(IIY) * RM;
    IDUM = (IA * IDUM + IC) % M;
    IIR[j] = IDUM;
    return randomValue;
}

void Static() {
    const int NV = 2;
    const double EPS = 1.0E-12f;
    const double RPI = 2.0f * asin(1.0f);

    static double Y[NV], DYDX[NV];
    static double RND[ny][nx], T[nPz], R[nPz], RK[nPz], DRK[nPz];
    static double IIR[97];

    std::string FNAME;
    double DZZ = 1.0f / float(nPz - 1);
    fill(&RU[0][0][0], 0.0f, nx * ny * nz);
    fill(&RV[0][0][0], 0.0f, nx * ny * nz);
    fill(&RW[0][0][0], 0.0f, nx * ny * nz);

    double SZ, SDZZDZ, SD2ZZDZ2;
#if !lShr
    #if lRem
    for (int k = 0; k < nPz; k++) {
        mkGrid<DIM::Z>(static_cast<float>(k) / static_cast<float>(nPz - 1), SZ, SDZZDZ, SD2ZZDZ2);
        getbackground(SZ, T[k], R[k]);
        RK[k] = kappa(SZ, R[k], T[k]);
        RK[k] *= 8.07f * THETA * REPR;
    }
    mkGrid<DIM::Z>(0, SZ, SDZZDZ, SD2ZZDZ2);
    DRK[0] = (-3.0f * RK[0] + 4.0f * RK[1] - RK[2]) * hz * SDZZDZ;
    for (int k = 1; k < nPz - 1; ++k) {
        mkGrid<DIM::Z>(static_cast<float>(k) / static_cast<float>(nPz - 1), SZ, SDZZDZ, SD2ZZDZ2);
        DRK[k] = (RK[k + 1] - RK[k - 1]) * hz * SDZZDZ;
    }
    mkGrid<DIM::Z>(1, SZ, SDZZDZ, SD2ZZDZ2);
    DRK[nPz - 1] = (3.0f * RK[nPz - 1] - 4.0f * RK[nPz - 2] + RK[nPz - 3]) * hz * SDZZDZ;

    if (myPE == 0) {
        int pos = fOut.find("  ") - 1;
        std::string FNAME = fOut.substr(0, pos) + ".strat";
        std::ofstream outFile(FNAME.c_str());
        for (int k = 0; k < nPz; ++k) {
            outFile << T[k] << " " << R[k] << " " << RK[k] / REPR << " " << DRK[k] / REPR << std::endl;
        }
        outFile.close();
    }

    TU = T[0];
    DZTB = T[nPz - 1] - 4.0f / 3.0f * T[nPz - 2] + 1.0f / 3.0f * T[nPz - 3];
    DZTU = T[0] - 4.0f / 3.0f * T[1] + 1.0f / 3.0f * T[2];
    #else
    T[0] = 1.0f, R[0] = 1.0f;
    Y[0] = 0.0f, Y[1] = 0.0f;

    double ZZ = 0.0f;

    for (int k = 1; k < nPz; ++k) {
        double HTRY = DZZ;
        derivs(ZZ, Y, NV, DYDX);
        double HDID, HNEXT;
        bsstep(Y, DYDX, NV, ZZ, HTRY, EPS, HDID, HNEXT);
        if (HDID != HTRY) exit_err("Static: Static structure error.");
        T[k] = 1.0f + Y[0];
        R[k] = exp(Y[1]) / (1.0f + Y[0]);
    }
    #endif
#else
    for (int k = 0; k < nPz; k++) {
        T[k] = 1.0f, R[k] = 1.0f;
    }
#endif
    // Divide among processors.
    if (myPEz == 0) {
        for (int z = 0; z < iz / 2; z++) {
            for (int x = 1; x < nx - ix + 1; x++) {
                for (int y = 1; y < ny - iy + 1; y++) {
                    TT[z][y][x] = T[0];
                    RO[z][y][x] = R[0];
                }
            }
#if lRem
            RKAPA[z] = RK[0];
            DKAPA[z] = DRK[0];
#endif
        }
        for (int z = iz / 2; z < nz; z++) {
            for (int x = 1; x < nx - ix + 1; x++) {
                for (int y = 1; y < ny - iy + 1; y++) {
                    TT[z][y][x] = T[z - iz / 2];
                    RO[z][y][x] = R[z - iz / 2];
                }
            }
#if lRem
            RKAPA[z] = RK[z - iz / 2];
            DKAPA[z] = DRK[z - iz / 2];
#endif
        }
    } else if (myPEz == nPEz - 1) {
        for (int z = 0; z < nRz + iz / 2; z++) {
            for (int y = 1; y < ny - iy + 1; y++) {
                for (int x = 1; x < nx - ix + 1; x++) {
                    TT[z][y][x] = T[nPz - nRz - iz / 2 + z];
                    RO[z][y][x] = R[nPz - nRz - iz / 2 + z];
                }
            }
#if lRem
            RKAPA[z] = RK[nPz - nRz - iz / 2 + z];
            DKAPA[z] = DRK[nPz - nRz - iz / 2 + z];
#endif
        }
        for (int z = nRz + iz / 2; z < nRz + iz; z++) {
            for (int y = 1; y < ny - iy + 1; y++) {
                for (int x = 1; x < nx - ix + 1; x++) {
                    TT[z][y][x] = T[nPz - 1];
                    RO[z][y][x] = R[nPz - 1];
                }
            }
#if lRem
            RKAPA[z] = RK[nPz - 1];
            DKAPA[z] = DRK[nPz - 1];
#endif
        }
    } else {
        int Start = myPEz * nRz - iz / 2;
        for (int z = 0; z < nz; z++) {
            for (int y = 1; y < ny - iy + 1; y++) {
                for (int x = 1; x < nx - ix + 1; x++) {
                    TT[z][y][x] = T[Start + z];
                    RO[z][y][x] = R[Start + z];
                }
            }
        }
#if lRem
        std::memcpy(RKAPA, RK + Start, sizeof(double) * nz);
        std::memcpy(DKAPA, DRK + Start, sizeof(double) * nz);
#endif
    }
#if lMag
    if (AMPB != 0.0) {
        double CLN = -4.0f * log(2.0f) / (BFH * BFH);
        for (int z = iz / 2 + 1; z <= nz - iz / 2; z++) {
            for (int y = 0; y < ny; y++) {
                for (int x = 0; x < nx; x++) {
                    BX[z][y][x] = AMPB * exp(CLN * square(ZEE[z] - BZP));
                }
            }
        }
        fill(&BY[0][0][0], 0.0f, nx * ny * nz);
        fill(&BZ[0][0][0], 0.0f, nx * ny * nz);
        for (int z = iz / 2 + 1; z <= nz - iz / 2; z++) {
            for (int y = 1; y < ny - iy + 1; y++) {
                for (int x = 1; x < nx - ix + 1; x++) {
                    RO[z][y][x] -= (BX[z][y][x] * BX[z][y][x] / TT[z][y][x]) * OBETA;
                }
            }
        }
    } else {
        tube();
    }
#endif
    if (AMPT != 0.0) {
#define ISW 1
#if ISW == 1
        int IDUM = -62659;
        for (int nnz = 0; nnz < nPEz; nnz++) {
            for (int k = iz / 2; k < nz - iz / 2; k++) {
                int iTag = k;
                for (int nny = 0; nny < nPEy; nny++) {
                    if (myPE == 0) {
                        for (int y = iy / 2; y < ny - iy / 2; y++) {
                            for (int x = ix / 2; x < nx - ix / 2; x++) {
                                RND[y][x] = 1.0f + AMPT * (RAN2(IDUM) - 0.5E00);
                            }
                        }
                        if ((nny + nPEy * nnz) != 0) {
                            MPI_Send(RND, nx * ny, MPI_DOUBLE, nny + nPEy * nnz, iTag, MPI_COMM_WORLD);
                        } else {
                            std::memcpy(&WW1[k][0][0], RND, sizeof(double) * nx * ny);
                        }
                    } else {
                        if (myPE == (nny + nPEy * nnz)) {
                            MPI_Status status;
                            MPI_Recv(&WW1[k][0][0], nx * ny, MPI_DOUBLE, 0, iTag, MPI_COMM_WORLD, &status);
                        }
                    }
                }
            }
        }

        for (int z = iz / 2; z < nz - iz / 2; z++) {
            for (int y = iy / 2; y < ny - iy / 2; y++) {
                for (int x = ix / 2; x < nx - ix / 2; x++) {
    #if !lShr
                    TT[z][y][x] = TT[z][y][x] * WW1[z][y][x];
    #else
                    RV[z][y][x] = RV[z][y][x] + WW1[z][y][x] - 1.0f;
    #endif
                }
            }
        }
#else
        double RPI = 2.0f * asin(1.0f);
        for (int z = iz / 2; z < nz - iz / 2; ++z) {
            for (int y = 1; y < ny - iy + 1; ++y) {
                for (int x = 1; x < nx - ix + 1; ++x) {
                    double RKY = WYY[y] / yMax * float(nPy) / (float(nPy) + 1.0f) * 2.0f * RPI * 8.0f;
                    double RKZ = ZEE[z] / zMax * 2.0f * RPI;
    #if !lShr
                    TT[z][y][x] = TT[z][y][x] + AMPT * sin(RKY) * sin(RKZ) / RO[z][y][x];
    #else
                    RV[z][y][x] = RV[z][y][x] + AMPT * sin(RKY) * sin(RKZ) / RO[z][y][x];
    #endif
                }
            }
        }
#endif
#undef ISW
    }
}

void tube() {
#if lMag
    fill(&BX[0][0][0], 0.0f, nx * ny * nz);
    fill(&BY[0][0][0], 0.0f, nx * ny * nz);
    fill(&BZ[0][0][0], 0.0f, nx * ny * nz);
#endif
    double OGAMMA = 1.0f / GAMMA;
    setup(ipar, par);
    static double XN[nTubes], ZN[nTubes], RRN[nTubes], expI[nTubes], HN[nTubes], QN[nTubes], BYN[nTubes], AN[nTubes];
    double XCENT, ZCENT, R_MAX, C_MT, LAMBDA, VZ0;
    int JCUT, xl, xr, yl, yr, zl, zr;
#if nTube == 0
    #if lMag
    fill(&BX[0][0][0], 0.0f, nx * ny * nz);
    fill(&BY[0][0][0], 0.0f, nx * ny * nz);
    fill(&BZ[0][0][0], 0.0f, nx * ny * nz);
    #endif
    fill(&RU[0][0][0], 0.0f, nx * ny * nz);
    fill(&RV[0][0][0], 0.0f, nx * ny * nz);
    fill(&RW[0][0][0], 0.0f, nx * ny * nz);
#elif nTube == 1
    XCENT = par[34];
    ZCENT = par[35];
    R_MAX = par[36];
    C_MT = par[37];
    A = par[38];
    C = exp(-R_MAX * R_MAX);
    double ROGAPR = pow(RO[2][1][2], GAMMA) / (RO[2][1][2] * TT[2][1][2]);
    JCUT = iy / 2 + 1;
    xl = 1, xr = nx - ix + 1;
    yl = 1, yr = ny - iy + 1;
    zl = iz / 2, zr = nz - iz / 2;
    for (int z = zl; z < zr; ++z) {
        for (int x = xl; x < xr; ++x) {
            double XNEW = EXX[x] - XCENT;
            double ZNEW = ZEE[z] - ZCENT;
            double RNEW = sqrt(XNEW * XNEW + ZNEW * ZNEW);

            if (RNEW < R_MAX) {
                BY[z][JCUT][x] = (exp(-RNEW * RNEW) - C) / (1.0f - C);
                double BPHI = BY[z][JCUT][x] * C_MT * A * pow(RNEW, 3) / (A * pow(RNEW, 3) + 1.0f);

                if (RNEW != 0.0) {
                    BX[z][JCUT][x] = BPHI * ZNEW / RNEW;
                    BZ[z][JCUT][x] = -BPHI * XNEW / RNEW;
                }

                double DPR_TOT;
                qromb(FF, RNEW, R_MAX, DPR_TOT);
                double DPR = OBETA / square(1.0f - C) * DPR_TOT - (square(BX[z][JCUT][x]) + square(BY[z][JCUT][x]) + square(BZ[z][JCUT][x])) * OBETA / 2.0f;
                double PRE = RO[z][JCUT][x] * TT[z][JCUT][x] + DPR;
                RO[z][JCUT][x] = pow(ROGAPR * PRE, OGAMMA);
                TT[z][JCUT][x] = PRE / RO[z][JCUT][x];
            }
        }
    }

    for (int z = zl; z < zr; ++z) {
        for (int y = yl; y < yr; ++y) {
            for (int x = xl; x < xr; ++x) {
                RO[z][y][x] = RO[z][JCUT][x];
                TT[z][y][x] = TT[z][JCUT][x];
                BX[z][y][x] = BX[z][JCUT][x];
                BY[z][y][x] = BY[z][JCUT][x];
                BZ[z][y][x] = BZ[z][JCUT][x];
            }
        }
    }
    fill(&RU[0][0][0], 0.0f, nx * ny * nz);
    fill(&RV[0][0][0], 0.0f, nx * ny * nz);
    fill(&RW[0][0][0], 0.0f, nx * ny * nz);
#elif nTube == 2
    XCENT = par[34];
    ZCENT = par[35];
    R_MAX = par[36];
    C_MT = par[37];
    A = par[38];
    C = exp(-R_MAX * R_MAX);
    double LAMBDA = par[51];
    double VZ0 = par[52];
    fill(&RW[0][0][0], 0.0f, nx * ny * nz);
    JCUT = iy / 2 + 1;
    xl = 1, xr = nx - ix + 1;
    yl = 1, yr = ny - iy + 1;
    zl = iz / 2, zr = nz - iz / 2;
    for (int z = zl; z < zr; ++z) {
        for (int x = xl; x < xr; ++x) {
            double XNEW = EXX[x] - XCENT;
            double ZNEW = ZEE[z] - ZCENT;
            double RNEW = sqrt(XNEW * XNEW + ZNEW * ZNEW);

            if (RNEW < R_MAX) {
                BY[z][JCUT][x] = (exp(-RNEW * RNEW) - C) / (1.0f - C);
                double BPHI = BY[z][JCUT][x] * C_MT * A * pow(RNEW, 3) / (A * pow(RNEW, 3) + 1.0f);

                if (RNEW != 0.0) {
                    BX[z][JCUT][x] = BPHI * ZNEW / RNEW;
                    BZ[z][JCUT][x] = -BPHI * XNEW / RNEW;
                }

                double DPR_TOT;
                qromb(FF, RNEW, R_MAX, DPR_TOT);
                double DPR = OBETA / square(1.0f - C) * DPR_TOT - (square(BX[z][JCUT][x]) + square(BY[z][JCUT][x]) + square(BZ[z][JCUT][x])) * OBETA / 2.0f;
                double PRE = RO[z][JCUT][x] * TT[z][JCUT][x] + DPR;
                TT[z][JCUT][x] = PRE / RO[z][JCUT][x];
            }
        }
    }

    for (int z = zl; z < zr; ++z) {
        for (int y = yl; y < yr; ++y) {
            for (int x = xl; x < xr; ++x) {
                TT[z][y][x] = TT[z][JCUT][x];
                BX[z][y][x] = BX[z][JCUT][x];
                BY[z][y][x] = BY[z][JCUT][x];
                BZ[z][y][x] = BZ[z][JCUT][x];
                RW[z][y][x] = RW[z][JCUT][x];
            }
        }
    }
    std::fill(&RU[0][0][0], &RU[0][0][0] + nx * ny * nz, 0.0f);
    std::fill(&RV[0][0][0], &RV[0][0][0] + nx * ny * nz, 0.0f);
    for (int z = zl; z < zr; ++z) {
        for (int y = yl; y < yr; ++y) {
            for (int x = xl; x < xr; ++x) {
                RW[z][y][x] *= sin(2 * PI * WYY[y] / LAMBDA - PI / 2);
            }
        }
    }
#elif nTube == 3
    XCENT = par[34];
    ZCENT = par[35];
    R_MAX = par[36];
    C_MT = par[37];
    A = par[38];
    C = exp(-R_MAX * R_MAX);
    LAMBDA = par[51];
    VZ0 = par[52];
    fill(&RW[0][0][0], 0.0f, nx * ny * nz);
    JCUT = iy / 2 + 1;
    xl = 1, xr = nx - ix + 1;
    yl = 1, yr = ny - iy + 1;
    zl = iz / 2, zr = nz - iz / 2;
    for (int z = zl; z < zr; ++z) {
        for (int x = xl; x < xr; ++x) {
            double XNEW = EXX[x] - XCENT;
            double ZNEW = ZEE[z] - ZCENT;
            double RNEW = sqrt(XNEW * XNEW + ZNEW * ZNEW);

            if (RNEW < R_MAX) {
                BY[z][JCUT][x] = (exp(-RNEW * RNEW) - C) / (1.0f - C);
                double BPHI = BY[z][JCUT][x] * C_MT * A * pow(RNEW, 3) / (A * pow(RNEW, 3) + 1.0f);

                if (RNEW != 0.0) {
                    BX[z][JCUT][x] = BPHI * ZNEW / RNEW;
                    BZ[z][JCUT][x] = -BPHI * XNEW / RNEW;
                }

                double DPR_TOT;
                qromb(FF, RNEW, R_MAX, DPR_TOT);
                double DPR = OBETA / square(1.0f - C) * DPR_TOT - (square(BX[z][JCUT][x]) + square(BY[z][JCUT][x]) + square(BZ[z][JCUT][x])) * OBETA / 2.0f;
                double PRE = RO[z][JCUT][x] * TT[z][JCUT][x] + DPR;
                TT[z][JCUT][x] = PRE / RO[z][JCUT][x];
            }
        }
    }

    for (int z = zl; z < zr; ++z) {
        for (int y = yl; y < yr; ++y) {
            for (int x = xl; x < xr; ++x) {
                TT[z][y][x] = TT[z][JCUT][x];
                BX[z][y][x] = BX[z][JCUT][x];
                BY[z][y][x] = BY[z][JCUT][x];
                BZ[z][y][x] = BZ[z][JCUT][x];
                RW[z][y][x] = RW[z][JCUT][x];
            }
        }
    }
    std::fill(&RU[0][0][0], &RU[0][0][0] + nx * ny * nz, 0.0f);
    std::fill(&RV[0][0][0], &RV[0][0][0] + nx * ny * nz, 0.0f);
    std::fill(&RW[0][0][0], &RW[0][0][0] + nx * ny * nz, 0.0f);
#elif nTube == 4
    ZCENT = par[35];
    R_MAX = par[36];
    C_MT = par[37];
    A = par[38];
    int JCUT = iy / 2 + 1;
    int xl = 1, xr = nx - ix + 1;
    int yl = 1, yr = ny - iy + 1;
    int zl = iz / 2, zr = nz - iz / 2;
    for (int z = zl; z < zr; ++z) {
        for (int x = xl; x < xr; ++x) {
            double ZNEW = ZEE[z] - ZCENT;
            double RNEW = sqrt(ZNEW * ZNEW);

            BY[z][JCUT][x] = (1.0f - tanh((RNEW - R_MAX) / A)) / 2;
            double PRE = RO[z][JCUT][x] * TT[z][JCUT][x] - BY[z][JCUT][x] * BY[z][JCUT][x] * OBETA / 2.0f;
            RO[z][JCUT][x] = PRE / TT[z][JCUT][x];
        }
    }
    for (int z = zl; z < zr; ++z) {
        for (int y = yl; y < yr; ++y) {
            for (int x = xl; x < xr; ++x) {
                BY[z][y][x] = BY[z][JCUT][x];
                RO[z][y][x] = RO[z][JCUT][x];
            }
        }
    }
    std::fill(&RU[0][0][0], &RU[0][0][0] + nx * ny * nz, 0.0f);
    std::fill(&RV[0][0][0], &RV[0][0][0] + nx * ny * nz, 0.0f);
#elif nTube == 5
    ZCENT = par[35];
    R_MAX = par[36];
    C_MT = par[37];
    A = par[38];
    double HH = par[55];
    bool ISTAG = 1;
    int N = 0;

    double DCOL = xMax / float(NCOL);
    double DROW = 2.0f * R_MAX / float(nRow);

    for (int k = 1; k <= nRow; k++) {
        for (int i = 1; i <= NCOL; i++) {
            if (ISTAG == 0) {
                XN[N - 1] = float(i) * DCOL - DCOL / 2.0f;
            } else {
                XN[N - 1] = float(i) * DCOL - float(k) * DCOL / 2.0f;
                if (XN[N - 1] < 0.0) {
                    XN[N - 1] = xMax + XN[N - 1];
                }
            }
            ZN[N - 1] = ZCENT - R_MAX - DROW / 2.0f + float(k) * DROW;
            N++;
        }
    }
    double DEN = -0.5f / square(HH);
    JCUT = iy / 2 + 1;
    xl = 1, xr = nx - ix + 1;
    yl = 1, yr = ny - iy + 1;
    zl = iz / 2, zr = nz - iz / 2;
    for (int z = zl; z < zr; z++) {
        for (int x = xl; x < xr; x++) {
            for (int N = 0; N < nTubes; N++) {
                double RR1 = fabs(EXX[x] - XN[N]);
                double RR2 = fabs(xMax - fabs(EXX[x] - XN[N]));
                if (RR1 <= xMax / 2.0f) {
                    BY[z][JCUT][x] += exp(DEN * RR1 * RR1) * exp(DEN * (ZEE[z] - ZN[N]) * (ZEE[z] - ZN[N]));
                }
                if (RR2 < xMax / 2.0f) {
                    BY[z][JCUT][x] += exp(DEN * RR2 * RR2) * exp(DEN * (ZEE[z] - ZN[N]) * (ZEE[z] - ZN[N]));
                }
            }
        }
    }
    exit_err("tube: Communication update to MPI needed");
        // TODO 4324 line in 3dmhdsub.f
#else
    exit_err("3dmhdtube: Invalid nTube number");
#endif
}

#endif
