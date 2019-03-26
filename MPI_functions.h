#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstdlib>

double **alloc_2d_int(int rows, int cols);

void padding (int *M_padded, int *N_padded, int M, int N, int sizei, int sizej);

void pack_buf (int Nx_p, int Ny_p, int Nx, int Ny, double **buf, double **masterbuf, int starti, int endi, int startj, int endj);
void unpack_buf (int Nx_p, int Ny_p, int Nx, int Ny, double **buf, double **masterbuf, int starti, int endi, int startj, int endj);




