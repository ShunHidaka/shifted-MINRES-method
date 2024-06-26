#include "../src/c/sminres_solver.h"
#include "functions.h"
char *FNAME = "./../data/sample-matrix.csr";

int M = 10;

int main(void){
  int i;
  /* read R-Symmetric or C-Hermitian Matrix 'A' */
  int N, DATASIZE;
  int *A_row, *A_col;
  double complex *A_ele;
  read_csr(FNAME, &N, &DATASIZE, &A_row, &A_col, &A_ele);
  /* generate shift points */
  double complex *sigma;
  sigma = (double complex *)calloc(M, sizeof(double complex));
  for(i=0; i<M; i++)
    sigma[i] = 0.01 * cexp( 2 * M_PI * I * (i + 0.5) / M );
  /* generate rhs-vector */
  double complex * b;
  b = (double complex *)calloc(N, sizeof(double complex));
  for(i=0; i<N; i++)
    b[i] = 1.0;

  int result;
  int itermax = 100000;
  double threshold = 1e-13;
  double complex *x;
  int status[3];
  x = (double complex *)calloc(M*N, sizeof(double complex));

  fprintf(stdout, "%s\n", FNAME);
  status[0]=1; status[1]=0; status[2]=20;
  result = solve_sminres(sigma, M, A_row,A_col,A_ele, x, b, N, itermax, threshold, status);

  if(result == 0)
    printf("fail\n");
  else
    printf("success\n");
  
  return 0;
}

/*
gcc sample1.c ../src/c/sminres_solver.c -lm -lgfortran -lblas -llapack
*/
