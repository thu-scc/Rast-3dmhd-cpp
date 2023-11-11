#ifndef _3DMHD_MPI_H_
#define _3DMHD_MPI_H_

#include "3dmhddata.h"

#include <mpi.h>
#include <vector>

void comm_mpi(data var);

void communicate() {
    comm_mpi(RU);
    comm_mpi(RV);
    comm_mpi(RW);
    comm_mpi(TT);
    comm_mpi(RO);
#if lMag
    comm_mpi(BX);
    comm_mpi(BY);
    comm_mpi(BZ);
#endif
}

MPI_Request reqs[4];
constexpr int BUFFER_Y_SIZE = nx * (iy / 2) * nz;
constexpr int BUFFER_Z_SIZE = nx * ny * (iz / 2);
using mpi_buffer = double[BUFFER_Y_SIZE];
// r : ny - iy  ->  ny - iy / 2
// l : 0        <-  iy / 2
mpi_buffer send_l, send_r, recv_l, recv_r;

void comm_mpi(data var) {
    int z, y, x;
    int xl = ix / 2, xr = nx - ix / 2;
    int yl = iy / 2, yr = ny - iy / 2;
    int zl = iz / 2, zr = nz - iz / 2;

    // Step 1.1: request MPI in Z direction
    if (nPEz > 1) {
        int iTagU = 100, iTagD = 200;
        if (myPEz != nPEz - 1) {
            MPI_Isend(&var[zr - zl][0][0], BUFFER_Z_SIZE, MPISIZE, myPE + nPEy, iTagU, MPI_COMM_WORLD, &reqs[0]);
            MPI_Irecv(&var[zr][0][0], BUFFER_Z_SIZE, MPISIZE, myPE + nPEy, iTagD, MPI_COMM_WORLD, &reqs[1]);
        } else {
            reqs[0] = MPI_REQUEST_NULL;
            reqs[1] = MPI_REQUEST_NULL;
        }

        if (myPEz != 0) {
            MPI_Irecv(&var[0][0][0], BUFFER_Z_SIZE, MPISIZE, myPE - nPEy, iTagU, MPI_COMM_WORLD, &reqs[2]);
            MPI_Isend(&var[zl][0][0], BUFFER_Z_SIZE, MPISIZE, myPE - nPEy, iTagD, MPI_COMM_WORLD, &reqs[3]);
        } else {
            reqs[2] = MPI_REQUEST_NULL;
            reqs[3] = MPI_REQUEST_NULL;
        }
    }

    // Step 1.2: prepare data in Y direction / swap data in Y direction
    if (nPEy > 1) {
        For(z, zl, zr) For(x, 0, nx) For(y, 0, yl) {
            send_r[(iy / 2) * (z * nx + x) + y] = var[z][ny - iy + y][x];
            send_l[(iy / 2) * (z * nx + x) + y] = var[z][iy / 2 + y][x];
        }
    } else {
        For(z, zl, zr) For(x, 0, nx) For(y, 0, yl) {
            var[z][y][x] = var[z][ny - iy + y][x];
            var[z][ny - iy / 2 + y][x] = var[z][iy / 2 + y][x];
        }
    }

    // Step 1.3: wait MPI in Z direction
    if (nPEz > 1) {
        MPI_Waitall(4, reqs, MPI_STATUS_IGNORE);
    }

    // Step 1.4: fix boundary in Y direction
    if (nPEy > 1) {
        For(z, 0, zl) For(x, 0, nx) For(y, 0, yl) {
            send_r[(iy / 2) * (z * nx + x) + y] = var[z][ny - iy + y][x];
            send_l[(iy / 2) * (z * nx + x) + y] = var[z][iy / 2 + y][x];
        }
        For(z, zr, nz) For(x, 0, nx) For(y, 0, yl) {
            send_r[(iy / 2) * (z * nx + x) + y] = var[z][ny - iy + y][x];
            send_l[(iy / 2) * (z * nx + x) + y] = var[z][iy / 2 + y][x];
        }
    } else {
        For(z, 0, zl) For(x, 0, nx) For(y, 0, yl) {
            var[z][y][x] = var[z][ny - iy + y][x];
            var[z][ny - iy / 2 + y][x] = var[z][iy / 2 + y][x];
        }
        For(z, zr, nz) For(x, 0, nx) For(y, 0, yl) {
            var[z][y][x] = var[z][ny - iy + y][x];
            var[z][ny - iy / 2 + y][x] = var[z][iy / 2 + y][x];
        }
    }

    // Step 2.1: request MPI in Y direction
    if (nPEy > 1) {
        int iTagR = 300, iTagL = 400;
        if (myPEy != nPEy - 1) {
            MPI_Isend(send_r, BUFFER_Y_SIZE, MPISIZE, myPE + 1, iTagR, MPI_COMM_WORLD, &reqs[0]);
            MPI_Irecv(recv_l, BUFFER_Y_SIZE, MPISIZE, myPE + 1, iTagL, MPI_COMM_WORLD, &reqs[1]);
        } else {
            MPI_Isend(send_r, BUFFER_Y_SIZE, MPISIZE, myPE - nPEy + 1, iTagR, MPI_COMM_WORLD, &reqs[0]);
            MPI_Irecv(recv_l, BUFFER_Y_SIZE, MPISIZE, myPE - nPEy + 1, iTagL, MPI_COMM_WORLD, &reqs[1]);
        }

        if (myPEy != 0) {
            MPI_Irecv(recv_r, BUFFER_Y_SIZE, MPISIZE, myPE - 1, iTagR, MPI_COMM_WORLD, &reqs[2]);
            MPI_Isend(send_l, BUFFER_Y_SIZE, MPISIZE, myPE - 1, iTagL, MPI_COMM_WORLD, &reqs[3]);
        } else {
            MPI_Irecv(recv_r, BUFFER_Y_SIZE, MPISIZE, myPE + nPEy - 1, iTagR, MPI_COMM_WORLD, &reqs[2]);
            MPI_Isend(send_l, BUFFER_Y_SIZE, MPISIZE, myPE + nPEy - 1, iTagL, MPI_COMM_WORLD, &reqs[3]);
        }
    }

    // Step 2.2: swap data in X direction
    For(z, 0, nz) For(y, yl, yr) For(x, 0, xl) {
        var[z][y][x] = var[z][y][xr - xl + x];
        var[z][y][xr + x] = var[z][y][xl + x];
    }

    // Step 2.3: wait MPI in Y direction
    if (nPEy > 1) {
        MPI_Waitall(4, reqs, MPI_STATUS_IGNORE);
    }

    // Step 3: copy data in Y direction, fix boundary in X direction
    if (nPEy > 1) {
        For(z, 0, nz) For(x, 0, nx) For(y, 0, yl) {
            var[z][ny - iy / 2 + y][x] = recv_l[(iy / 2) * (z * nx + x) + y];
            var[z][y][x] = recv_r[(iy / 2) * (z * nx + x) + y];
        }
    }
    For(z, 0, nz) For(y, 0, yl) For(x, 0, xl) {
        var[z][y][x] = var[z][y][xr - xl + x];
        var[z][y][xr + x] = var[z][y][xl + x];
    }
    For(z, 0, nz) For(y, yr, ny) For(x, 0, xl) {
        var[z][y][x] = var[z][y][xr - xl + x];
        var[z][y][xr + x] = var[z][y][xl + x];
    }
}

#endif
