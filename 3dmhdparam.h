#pragma once
#include <string>
#include <sys/stat.h>
const std::string outputDir = "output/";
const std::string fInp = outputDir + "test";
const std::string fOut = outputDir + "test";

/* `lRot`: Rotation */
#define lRot false
/* `lMag`: Magentic fields */
#define lMag false
/* `lPot`: Ensure divB=0 */
#define lPot false
/* `lRem`: Rempel relax */
#define lRem false
/* `lShr`: Shear instability */
#define lShr false

const int nPEy = 4;
/* The minimum # of processors in the z direction is 2 at the moment! */
const int nPEz = 6;
const int nPE = nPEy * nPEz;

/* `nPx`, `nPy` and `nPz` are the number of grid points in each direction. */
const int nPx = 256;
const int nPy = 256;
const int nPz = 256;

/* `nRz` is the number of vertical grid points resident on each processor. */
const int nRy = nPy / nPEy;
const int nRz = nPz / nPEz;

/* `nCore` is the number of grid points resident on each processor. */
const int nCore = nPx * nRy * nRz;

// 必须是偶数，否则计算会出现 bug
/* `ix` and `iy` must be at least two for second-order horizontal derivatives. */
const int ix = 2;
const int iy = 2;
/* `iz=2` for second-order vertical derivatives. */
const int iz = 2;
/* `IPAD` pads the memory location to increase perfomance. */
const int IPAD = 0;

const int ny = nRy + iy; // R
const int nz = nRz + iz; // R
const int nx = nPx + ix; // P

/**
 * The function used for the stretching of the grid is determined by the parameter `nGrid`.
 * Putting `nGrid=1` or `nGrid=2` will make the code going through different pieces of the routine `mkGrid`
 *   which is called by `xJacobi`, `yJacobi` and `zJacobi`.
 * To add a new grid definition just edit `mkGrid`.
 *
 * - `nGrid=0`: Uniform grid.
 *
 * - `nGrid=1`: Original arctan grid function.
 *    `XX1`, `XX2`, `YY1`, `YY2`, and `ZZ1`, `ZZ2` are the arctan parameters for the grid stretch
 *     (set to very small values, -1.0E-09, 1.0E-09 for no stretching.
 *
 * - `nGrid=2`: A function containing two ArcTan is used for the stretching.
 *    It produces a dip in grid spacing.
 *     `XA`: determines the position of the dip. Takes a value between 0. and 1.
 *     `XB`: determines the half width of the dip. 0.2 means a width of 0.4.
 *     `XC`: determines the sharpness of the walls. Set to 0 to have a regular grid.
 *     `XD`: determines the depth of the dip in terms of DDX(center of dip)/DDX(0). must be between 1.E-9 and 1-1.E-9.
 */
#define nGrid 0

/* for `nGrid=1`, see also `nGrid` and `mkGrid` */
const double XX1 = 0, XX2 = 0;
const double YY1 = 0, YY2 = 0;
const double ZZ1 = 0, ZZ2 = 0;

/* for `nGrid=2`, see also `nGrid` and `mkGrid` */
const double XA = 0, XB = 0, XC = 0, XD = 0;
const double YA = 0, YB = 0, YC = 0, YD = 0;
const double ZA = 0, ZB = 0, ZC = 0, ZD = 0;

/**
 * `ixCon` and `iyCon` specify type of horizontal boundary condition.
 * - =0 periodic.
 * NOTE: Horizontal periodicity hard-wired into border values in subroutine `step`.
 * NOTE: Upper and lower boundaries hard-wired to be stress free and impenetrable in subroutine `BCON`.
 **/
#define ixCon 0
#define iyCon 0
/**
 * `izCon` specifies lower boundary temperature condition.
 * - =0 constant temperature.
 * - =1 constant thermal flux.
 **/
#define izCon 1
/**
 * `iTCon` specifies plume source and/or upper boundary temperature condition.
 * - =0 fixed temperature.
 * - =1 fixed flux.
 * - =2 embedded heat loss.
 * - =3 thermal.
 * - =4 nonuniform boundary.
 **/
#define iTCon 4
/**
 * `iBCon` specifies the upper and lower boundary magnetic field conditions.
 * - =0 Bz=0
 **/
#define iBCon 0

// const double SFF = 0.315f;
// const double SFF = 0.4f;
const double SFF = 0.8f;
/**
 * `iD` specifies the vertical variation in density diffusion.
 * - =0 no diffusion.
 * - =1 goes as 1/rhobar.
 * - =2 exponentially confined to upper and lower boundaries.
 */
#define iD 0
const double DH = 0, DV = 0;
const int iWord = 8;

const int nCol = 0, nRow = 0, nTubes = nCol * nRow;

#define nTube 0