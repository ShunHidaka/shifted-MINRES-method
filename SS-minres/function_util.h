#include "definition.h"

void solve_sminres(double complex *sigma,
                   int *A_row, int *A_col, double complex *A_ele,
                   double complex *rhs,
                   double complex *x, int M, int N,
                   int itermax, double threshold, int *status);

void SpMV(const int *A_row, const int *A_col, const double complex *A_ele,
          const double complex *x, double complex *b, int N);

FILE *fopen_mtx(const char *fname, const char *mode, int *row_size, int *col_size, int *ele_size);


void read_csr(const char *fname, int *N, int *DATASIZE,
              int **row_ptr, int **col_ind, double complex **element);
