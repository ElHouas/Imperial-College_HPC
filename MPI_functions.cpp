#include "MPI_functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdexcept>
#include <iostream>
#include <cstdlib>

double **alloc_2d_int(int rows, int cols) {
    double *data = (double *)malloc(rows*cols*sizeof(double));
    double **array= (double **)malloc(rows*sizeof(double*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[i*cols]);
    
    for (int i=0; i < rows; i++)
        {
            for (int j=0; j < cols; j++)
            {
               array[i][j]=0;
            }
        } 
    
    return array;
}


/* This function sets M_padded (or N_padded) to a size divisible by p */

void padding (int *M_padded, int *N_padded, int M, int N, int sizei, int sizej)
{
  if (M%sizei != 0)
  {
    *M_padded = M + sizei - M % sizei;
  }
  else
  {
    *M_padded = M;
  }
  if (N%sizej != 0)
  {
    *N_padded = N + sizej - N % sizej;
  }
  else
  {
    *N_padded = N;
  }
}

 /* Set buf to appropriate values of masterbuf for each process */
void pack_buf (int Nx_p, int Ny_p, int Nx, int Ny, double **buf, double **masterbuf, int starti, int endi, int startj, int endj)
{
  int i, j;
  int x = 0;
  int y = 0;
  for (unsigned i = starti; i < endi; i++)
  {
    y = 0;
    for (unsigned j = startj; j < endj; j++)
    {
        buf[x][y] = masterbuf[i][j];
      y++;
    }
  x++;
  }

}


  /* Set appropriate values of masterbuf to values of buf for each process */
void unpack_buf (int Nx_p, int Ny_p, int Nx, int Ny, double **buf, double **masterbuf, int starti, int endi, int startj, int endj)
{
  int i, j;
  int x = 0;
  int y = 0;
  for (i = starti; i < endi; i++)
  {
    y = 0;
    for (j = startj; j < endj; j++)
    {
      masterbuf[i][j] = buf[x][y];
      y++;
    }
  x++;
  }
}