#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include "function_util.h"
#include "function_blas.h"

char *FNAME;
int MAX_ITR = INT_MAX;
int ZERO = 0;
int ONE = 1;
double dTMP = 0;
double complex cTMP = 0.0 + 0.0I;

void set_fname(int c, char *v[])
{
  FNAME = "./../../MATRIX/PPE3594_A.csr";
}

void set_shifts(int *M, double complex **sigma)
{
  *M = 4;
  *sigma = (double complex *)calloc(*M, sizeof(double complex));
  (*sigma)[0] = 0.0 + 0.00I;
  (*sigma)[1] = 0.0 + 0.01I;
  (*sigma)[2] = 0.0 + 0.10I;
  (*sigma)[3] = 0.0 + 1.00I;
}

void set_rhsVector(int N, double complex **b, double *norm)
{
  *b = (double complex *)calloc(N, sizeof(double complex));
  for(int i=0; i<N; i++)
    (*b)[i] = 1.0 + 0.0I;
  *norm = dznrm2_(&N, *b, &ONE);
}

void SpMV(const int *A_row, const int *A_col, const double complex *A_ele,
	  const double complex *x, double complex *b, int N)
{
  int i, j;
  double complex tmp;
#pragma omp parallel for private(j, tmp)
  for(i=0; i<N; i++){
    tmp = 0.0+0.0I;
    for(j=A_row[i]; j<A_row[i+1]; j++)
      tmp += A_ele[j]*x[A_col[j]];
    b[i] = tmp;
  }
}


FILE *fopen_mtx(const char *fname, const char *mode, int *row_size, int *col_size, int *ele_size)
{
  FILE *fp = NULL;
  fp = fopen(fname, mode);
  if(fp == NULL){
    fprintf(stderr, "Can not open file : %s\n", fname);
    exit(1);
  }
  int chr;
  while( (chr = fgetc(fp)) != EOF && (chr=='%' || chr=='#') )
    while( (chr = fgetc(fp)) != EOF )
      if(chr == '\n') break;
  fseek(fp, -sizeof(char), SEEK_CUR);
  int num, tmp1, tmp2, tmp3;
  num = fscanf(fp, "%d %d %d", &tmp1, &tmp2, &tmp3);
  if(num != 3){ fprintf(stderr, "fopen err\n"); exit(1);}
  if(col_size != NULL) *row_size = tmp1;
  if(row_size != NULL) *col_size = tmp2;
  if(ele_size != NULL) *ele_size = tmp3;
  return fp;
}
void read_csr(const char *fname, int *N, int *DATASIZE,
              int **row_ptr, int **col_ind, double complex **element)
{
  FILE *fp;
  int rsize, csize, esize;
  fp = fopen_mtx(fname, "r", &rsize, &csize, &esize);
  *N = rsize-1; *DATASIZE = esize;

  *row_ptr = (int *)calloc(rsize, sizeof(int));
  *col_ind = (int *)calloc(csize, sizeof(int));
  *element = (double complex *)calloc(esize, sizeof(double complex));

  int row, col, i=0;
  double real, imag=0;
  while( fscanf(fp, "%d %d %lf %lf", &row, &col, &real, &imag) != EOF ){
    if(i < rsize) (*row_ptr)[i] = row;
    if(i < csize) (*col_ind)[i] = col;
    if(i < esize) (*element)[i] = real + imag*I;
    i = i + 1;
  }
  fclose(fp);
}
