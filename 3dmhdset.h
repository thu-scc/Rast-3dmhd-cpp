#ifndef _3DMHD_SET_H_
#define _3DMHD_SET_H_

#include <cmath>
#include <string>

void setup(int ipar[32], double par[64]) {
    for (int i = 0; i < 32; i++) ipar[i] = 0;
    for (int i = 0; i < 64; i++) par[i] = 0;
    ipar[0] = 1;
    ipar[1] = 1;
    ipar[2] = 10;
    ipar[3] = 1;
    par[0] = 75;
    par[1] = 0.2f;
    par[2] = 5;
    par[3] = 10;
    par[11] = 0.1f;
    par[16] = 6.0f;
    par[17] = 6.0f;
    par[18] = 2.0f;
    par[19] = 1.0e-03f;
    par[53] = 5.0f / 3.0f;
}

/**
 * defination of `ipar`
 * --------------------------------------------------------------
 *   00  nCase  (case number)
 *   01  nCaseP (previous case number if continuation run)
 *   02  nTotal (total number of time steps required)
 *   03  nStep0 (number of steps between field '.dat0','.par',and '.lis' dumps)
 *   04  nStart (0=new start, else= dump number from previous solution)
 *   05  ixCon  (adjustable in parameter statement, x boundary condition)
 *   06  iyCon  (adjustable in parameter statement, y boundary condition)
 *   07  izCon  (adjustable in parameter statement, lower boundary condition)
 *   08  iTCon  (adjustable in parameter statement, upper temperature condition)
 *   09  iBCon  (adjustable in parameter statement, magnetic field condition)
 *   10  nPx    (adjustable in parameter statement, number x grid points)
 *   11  nPy    (adjustable in parameter statement, number x grid points)
 *   12  nPz    (adjustable in parameter statement, number z grid points)
 *   13  nPE    (adjustable in parameter statement, total number of processors)
 *   14  nIt    (nonadjustable, time step number)
 *   15  nTube  (adjustable in parameter statement, tube's model number (initial condition))
 *   16  nGrid  (adjustable in parameter statement, grid number (type of grid))
 *   17  nCol   (adjustable in parameter statement, number of columns of tubes)
 *   18  nRow   (adjustable in parameter statement, number of rows of tubes)
 *   19  iD     (adjustable in parameter statement, density diffusion model)
 *   20  nPEy   (adjustable in parameter statement, processors in y direction)
 * --------------------------------------------------------------
 **/

/**
 * defination of `par`
 * --------------------------------------------------------------
 *   00  RE    (soundspeed Reynolds number)
 *   01  PR    (pseudo Prandtl number)
 *   02  THETA (temperature gradient, THETA=ETA for lRem (M. Rempel))
 *   03  GRAV  (nondimensional gravitational constant, = (M+1)*THETA; M=1.5)
 *   04  RY    (soundspeed Rossby number)
 *   05  ANG   (rotation axis angle from vertical)
 *   06  RM    (magnetic Reynolds number)
 *   07  BETA  (plasma beta)
 *   08  PZP   (depth of transition to stable background)
 *   09  SIGMA (width of transition region)
 *   10  POLYS (polytropic index of lower stable layer)
 *   11  TP    (plume temperature or flux, or boundary nonuniformity)
 *   12  TC    (cooling time scale for perturbation onset; embedded onset time)
 *   13  XP    (plume source x position)
 *   14  YP    (plume source y position)
 *   15  ZP    (plume source z position)
 *   16  XMAX  (maximum x dimension)
 *   17  YMAX  (maximum y dimension)
 *   18  ZMAX  (maximum z dimension)
 *   19  AMPT  (amplitude of new start random temperature fluctuations)
 *   20  AMPB  (amplitude of new start magnetic field layer)
 *   21  BFH   (full-width at half-maximum of new start magnetic field layer)
 *   22  BZP   (depth of maximum of new start magnetic field layer)
 *   23  XX1   (adjustable in parameter statement, grid stretching parameter)
 *   24  XX2   (adjustable in parameter statement, grid stretching parameter)
 *   25  YY1   (adjustable in parameter statement, grid stretching parameter)
 *   26  YY2   (adjustable in parameter statement, grid stretching parameter)
 *   27  ZZ1   (adjustable in parameter statement, grid stretching parameter)
 *   28  ZZ2   (adjustable in parameter statement, grid stretching parameter)
 *   29  SFF   (adjustable in parameter statement, time stepping safety factor)
 *   30  DT    (nonadjustable, time step)
 *   31  TIMT  (nonadjustable, total simulation time)
 *   32  TIMC  (nonadjustable, current simulation time)
 *   33  XCENT (initial x position of the center of the magnetic tube)
 *   34  ZCENT (initial z position of the center of the magnetic tube)
 *   35  RMAX  (radius at which the magnetic field is set to zero)
 *   36  CMT   (determine the pitch angle where B_phi has its maximum)
 *   37  A     (determine the profile of B_phi)
 *   38  XA    (adjustable in parameter statement, grid stretching parameter)
 *   39  XB    (adjustable in parameter statement, grid stretching parameter)
 *   40  XC    (adjustable in parameter statement, grid stretching parameter)
 *   41  XD    (adjustable in parameter statement, grid stretching parameter)
 *   42  YA    (adjustable in parameter statement, grid stretching parameter)
 *   43  YB    (adjustable in parameter statement, grid stretching parameter)
 *   44  YC    (adjustable in parameter statement, grid stretching parameter)
 *   45  YD    (adjustable in parameter statement, grid stretching parameter)
 *   46  ZA    (adjustable in parameter statement, grid stretching parameter)
 *   47  ZB    (adjustable in parameter statement, grid stretching parameter)
 *   48  ZC    (adjustable in parameter statement, grid stretching parameter)
 *   49  ZD    (adjustable in parameter statement, grid stretching parameter)
 *   50  LAMBDA(wavelength of the initial pertubation of the 3d tube)
 *   51  VZ0   (amplitude  of the initial pertubation of the 3d tube)
 *   52  QFH   (temporal full-width at half-maximum of embedded heat loss)
 *   53  GAMMA (ratio of specific heats: Cp/Cv)
 *   54  HH    (tube widths for layer of tubes, spatial width for plume)
 *   55  DH    (adjustable in parameter statement, horizontal density diffusion)
 *   56  DV    (adjustable in parameter statement, vertical density diffusion)
 *   57  TSTART(time for heat flux relaxation start)
 *   58  TOFF  (time scale for heat flux relaxation termination)
 *   59  RLAX  (amplitude of heat flux relaxation)
 * --------------------------------------------------------------
 **/

#endif
