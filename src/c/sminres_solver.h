#include "definition.h"

// CSR形式の行列Aおよびシフトsigmaを渡して近似解を返す
int solve_sminres(const double complex *sigma, const int M,
                  const int *A_row, const int *A_col, const double complex *A_ele,
                  double complex *x, double complex *rhs, int N,
                  int itermax, double threshold, int *status);

// solve_sminres()内で用いる行列ベクトル積関数
void sample_SpMV(const int *A_row, const int *A_col, const double complex *A_ele,
                 const double complex *x, double complex *b, int N);
