#ifndef _3DMHD_UTILS_H_
#define _3DMHD_UTILS_H_

#include "3dmhdparam.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

const double PI = 2.0f * asin(1.0f);
constexpr double DOUBLE_INF = 1e99;

#define For(i_, l_, r_) for (i_ = l_; i_ < r_; i_++)
#define FOR3D           For(z, zl, zr) For(y, yl, yr) For(x, xl, xr)
#define SingleFOR3D \
    int z, y, x;    \
    FOR3D
#define FORRANGE      zl, zr, yl, yr, xl, xr
#define FORRANGEPARAM int zl, int zr, int yl, int yr, int xl, int xr
#define FORRANGE_EXT  zl - 1, zr  + 1, yl - 1, yr + 1, xl - 1, xr + 1
#define FORRANGE_EXT_XY zl, zr, yl - 1, yr + 1, xl - 1, xr + 1
#define FORRANGE_EXT_X zl, zr, yl, yr, xl - 1, xr + 1
#define FORRANGE_EXT_Y zl, zr, yl - 1, yr + 1, xl, xr

using std::max;
using std::min;

void standard(char* str) {
    bool flag = false;
    for (int i = 0; i < 79; ++i) {
        if (str[i] == '@') flag = true;
        if (flag) str[i] = ' ';
    }
    str[79] = '\n';
}

void exit_err(std::string info) {
    std::cerr << info << std::endl;
    MPI_Finalize();
    exit(-1);
}

void fill(double* src, double val, int n) {
    std::fill(src, src + n, val);
}

void fill(data src, double val, FORRANGEPARAM) {
    SingleFOR3D src[z][y][x] = val;
}

void copy(data tar, data src, FORRANGEPARAM) {
    SingleFOR3D tar[z][y][x] = src[z][y][x];
}

void scale_add(data tar, data src1, data srcd1, double scale, FORRANGEPARAM) {
    SingleFOR3D tar[z][y][x] = src1[z][y][x] + scale * srcd1[z][y][x];
}

double minVal(double* src, int l, int r) {
    double res = DOUBLE_INF;
    for (int i = l; i < r; i++)
        res = min(res, src[i]);
    return res;
}

double maxVal(double* src, int l, int r) {
    double res = -DOUBLE_INF;
    for (int i = l; i < r; i++)
        res = max(res, src[i]);
    return res;
}

double minVal(data src, FORRANGEPARAM) {
    double res = DOUBLE_INF;
    SingleFOR3D res = min(res, src[z][y][x]);
    return res;
}

void minLoc(data src, int a[3], FORRANGEPARAM) {
    double res = DOUBLE_INF;
    SingleFOR3D {
        if (res > src[z][y][x]) {
            a[0] = x, a[1] = y, a[2] = z;
            res = src[z][y][x];
        }
    }
}

double maxVal(data src, FORRANGEPARAM) {
    double res = -DOUBLE_INF;
    SingleFOR3D res = max(res, src[z][y][x]);
    return res;
}

double buffer[nCore];

void load(std::ifstream& in, data src, FORRANGEPARAM) {
    int size = (xr - xl) * (yr - yl) * (zr - zl);
    in.read(reinterpret_cast<char*>(buffer), size * sizeof(double));
    SingleFOR3D src[z][y][x] = buffer[(x - xl) + (xr - xl) * (y - yl) + (xr - xl) * (yr - yl) * (z - zl)];
}

void save(std::ofstream& out, data src, FORRANGEPARAM) {
    int size = (xr - xl) * (yr - yl) * (zr - zl);
    assert(size == nCore);
    SingleFOR3D buffer[(x - xl) + (xr - xl) * (y - yl) + (xr - xl) * (yr - yl) * (z - zl)] = src[z][y][x];
    out.write(reinterpret_cast<char*>(buffer), size * sizeof(double));
}

enum OPT { SET = 0,
           ADD = 1,
           SUB = 2 };
enum DIM { Z = 0,
           Y = 1,
           X = 2 };

#define ci constexpr inline
#define defineDimSwitch(type, u, x, y, z)         \
    template <DIM> ci type u();                   \
    template <> ci type u<DIM::X>() { return x; } \
    template <> ci type u<DIM::Y>() { return y; } \
    template <> ci type u<DIM::Z>() { return z; }

template <OPT> constexpr inline double& oper(double& a, double b);
template <> constexpr inline double& oper<OPT::SET>(double& a, double b) { return a = b; }
template <> constexpr inline double& oper<OPT::ADD>(double& a, double b) { return a += b; }
template <> constexpr inline double& oper<OPT::SUB>(double& a, double b) { return a -= b; }

template <DIM> constexpr inline const double& extend(const double* a, int z, int y, int x);
template <> constexpr inline const double& extend<DIM::X>(const double* a, int z, int y, int x) { return a[x]; }
template <> constexpr inline const double& extend<DIM::Y>(const double* a, int z, int y, int x) { return a[y]; }
template <> constexpr inline const double& extend<DIM::Z>(const double* a, int z, int y, int x) { return a[z]; }

// d[f]/dx (x) = (f(x+h) - f(x-h)) / 2h
template <DIM> constexpr inline double dFdu(const data F, int z, int y, int x);
template <> constexpr inline double dFdu<DIM::X>(const data F, int z, int y, int x) { return F[z][y][x + 1] - F[z][y][x - 1]; }
template <> constexpr inline double dFdu<DIM::Y>(const data F, int z, int y, int x) { return F[z][y + 1][x] - F[z][y - 1][x]; }
template <> constexpr inline double dFdu<DIM::Z>(const data F, int z, int y, int x) { return F[z + 1][y][x] - F[z - 1][y][x]; }

template <DIM> constexpr inline double dFGdu(const data F, const data G, int z, int y, int x);
template <> constexpr inline double dFGdu<DIM::X>(const data F, const data G, int z, int y, int x) { return F[z][y][x + 1] * G[z][y][x + 1] - F[z][y][x - 1] * G[z][y][x - 1]; }
template <> constexpr inline double dFGdu<DIM::Y>(const data F, const data G, int z, int y, int x) { return F[z][y + 1][x] * G[z][y + 1][x] - F[z][y - 1][x] * G[z][y - 1][x]; }
template <> constexpr inline double dFGdu<DIM::Z>(const data F, const data G, int z, int y, int x) { return F[z + 1][y][x] * G[z + 1][y][x] - F[z - 1][y][x] * G[z - 1][y][x]; }

template <DIM> constexpr inline const double& duudu(int z, int y, int x);
template <> constexpr inline const double& duudu<DIM::X>(int z, int y, int x) { return dxxdx[x]; }
template <> constexpr inline const double& duudu<DIM::Y>(int z, int y, int x) { return dyydy[y]; }
template <> constexpr inline const double& duudu<DIM::Z>(int z, int y, int x) { return dzzdz[z]; }

defineDimSwitch(const double&, hu, hx, hy, hz);

// d^2[f]/dx^2 (x) = (f(x+h) - 2f(x) + f(x-h)) / h^2
template <DIM>
constexpr inline double d2Fdu2(const data a, int z, int y, int x);
template <> constexpr inline double d2Fdu2<DIM::X>(const data a, int z, int y, int x) { return a[z][y][x + 1] - 2 * a[z][y][x] + a[z][y][x - 1]; }
template <> constexpr inline double d2Fdu2<DIM::Y>(const data a, int z, int y, int x) { return a[z][y + 1][x] - 2 * a[z][y][x] + a[z][y - 1][x]; }
template <> constexpr inline double d2Fdu2<DIM::Z>(const data a, int z, int y, int x) { return a[z + 1][y][x] - 2 * a[z][y][x] + a[z - 1][y][x]; }

template <DIM> constexpr inline const double& d2uudu2(int z, int y, int x);
template <> constexpr inline const double& d2uudu2<DIM::X>(int z, int y, int x) { return d2xxdx2[x]; }
template <> constexpr inline const double& d2uudu2<DIM::Y>(int z, int y, int x) { return d2yydy2[y]; }
template <> constexpr inline const double& d2uudu2<DIM::Z>(int z, int y, int x) { return d2zzdz2[z]; }

defineDimSwitch(const double&, h2u, h2x, h2y, h2z);

template <OPT opt>
void func_base(FORRANGEPARAM, data tar, const data src, double scale = 1) {
    SingleFOR3D oper<opt>(tar[z][y][x], src[z][y][x] * scale);
}

template <OPT opt, DIM dim> void func_base_extend(FORRANGEPARAM, data tar, const double* src, double scale = 1) {
    SingleFOR3D oper<opt>(tar[z][y][x], extend<dim>(src, z, y, x) * scale);
}

template <OPT opt> void func_mul2(FORRANGEPARAM, data tar, const data src1, const data src2, double scale = 1) {
    SingleFOR3D oper<opt>(tar[z][y][x], src1[z][y][x] * src2[z][y][x] * scale);
}

template <OPT opt> void func_mul3(FORRANGEPARAM, data tar, const data src1, const data src2, const data src3, double scale = 1) {
    SingleFOR3D oper<opt>(tar[z][y][x], src1[z][y][x] * src2[z][y][x] * src3[z][y][x] * scale);
}

template <OPT opt, DIM dim> void func_mul_extend(FORRANGEPARAM, data tar, const data src1, const double* src2, double scale = 1) {
    SingleFOR3D oper<opt>(tar[z][y][x], src1[z][y][x] * extend<dim>(src2, z, y, x) * scale);
}

void func_inv(FORRANGEPARAM, data tar) {
    SingleFOR3D tar[z][y][x] = 1.0f / tar[z][y][x];
}

// TODO SIMD

template <OPT opt, DIM dim> void func_d1(FORRANGEPARAM, data tar, const data src, double scale = 1) {
    SingleFOR3D oper<opt>(tar[z][y][x],
                          dFdu<dim>(src, z, y, x) * duudu<dim>(z, y, x) * hu<dim>() *
                              scale);
}

template <OPT opt, DIM dim> void func_d1_mul(FORRANGEPARAM, data tar, const data src1, const data src2, double scale = 1) {
    SingleFOR3D oper<opt>(tar[z][y][x],
                          dFdu<dim>(src1, z, y, x) * duudu<dim>(z, y, x) * hu<dim>() *
                              src2[z][y][x] * scale);
}

template <OPT opt, DIM dim> void func_d1_mul_array_z(FORRANGEPARAM, data tar, const data src1, const data src2, const double* array_z, double scale = 1) {
    SingleFOR3D oper<opt>(tar[z][y][x],
                          dFdu<dim>(src1, z, y, x) * duudu<dim>(z, y, x) * hu<dim>() *
                              src2[z][y][x] * array_z[z] * scale);
}

constexpr inline double square(double x) { return x * x; }

template <OPT opt, DIM dim> void func_d2(FORRANGEPARAM, data tar, const data src, double scale = 1) {
    SingleFOR3D oper<opt>(tar[z][y][x],
                          (dFdu<dim>(src, z, y, x) * d2uudu2<dim>(z, y, x) * hu<dim>() +
                           d2Fdu2<dim>(src, z, y, x) * square(duudu<dim>(z, y, x)) * h2u<dim>()) *
                              scale);
}

template <OPT opt, DIM dim> void func_d2_mul(FORRANGEPARAM, data tar, const data src1, const data src2, double scale = 1) {
    SingleFOR3D oper<opt>(tar[z][y][x],
                          (dFdu<dim>(src1, z, y, x) * d2uudu2<dim>(z, y, x) * hu<dim>() +
                           d2Fdu2<dim>(src1, z, y, x) * square(duudu<dim>(z, y, x)) * h2u<dim>()) *
                              src2[z][y][x] * scale);
}

template <OPT opt, DIM dim> void func_d2_mul_array_z(FORRANGEPARAM, data tar, const data src1, const data src2, const double* array_z, double scale = 1) {
    SingleFOR3D oper<opt>(tar[z][y][x],
                          (dFdu<dim>(src1, z, y, x) * d2uudu2<dim>(z, y, x) * hu<dim>() +
                           d2Fdu2<dim>(src1, z, y, x) * square(duudu<dim>(z, y, x)) * h2u<dim>()) *
                              src2[z][y][x] * array_z[z] * scale);
}

template <OPT opt, DIM dim> void func_mul_d1(FORRANGEPARAM, data tar, const data src1, const data src2) {
    SingleFOR3D oper<opt>(tar[z][y][x],
                          dFGdu<dim>(src1, src2, z, y, x) * duudu<dim>(z, y, x) * hu<dim>());
}

#endif
