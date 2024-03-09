#include "definition.h"

void sminres_init();
void sminres_update();
void sminres_finalize();

int solve_sminres(const double complex *sigma, const int M,
                  const int *A_row, const int *A_col, const double complex *A_ele,
                  double complex *x, double complex *rhs, int N,
                  int itermax, double threshold, int *status);
void sample_SpMV(const int *A_row, const int *A_col, const double complex *A_ele,
                 const double complex *x, double complex *b, int N);

FILE *fopen_mtx(const char *fname, const char *mode, int *row_size, int *col_size, int *ele_size);
void read_csr(const char *fname, int *N, int *DATASIZE,
              int **row_ptr, int **col_ind, double complex **element);
