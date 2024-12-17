#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

FILE *fopen_mtx(const char *fname, const char *mode, int *row_size, int *col_size, int *ele_size);
void read_csr(const char *fname, int *N, int *DATASIZE, int **row_ptr, int **col_ind, double complex **element);

void SpMV(const int *A_row, const int *A_col, const double complex *A_ele,
          const double complex *x, double complex *b, int N);
