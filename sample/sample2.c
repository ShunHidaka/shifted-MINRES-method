#include "../src/c/sminres_function.h"
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
  /* */
  int itermax = 100000;
  double threshold = 1e-13;
  double complex *x, *v, *Av;
  int status;
  x  = (double complex *)calloc(M*N, sizeof(double complex));
  v  = (double complex *)calloc(N,   sizeof(double complex));
  Av = (double complex *)calloc(N,   sizeof(double complex));

  sminres_init(N, b, M, sigma, v, x, threshold);
  for(i=0; i<itermax; i++){
    SpMV(A_row,A_col,A_ele, v, Av, N);
    sminres_update(i, v, Av, x, &status);
    if(status == 1){
      printf("Converged %d\n", i);
      break;
    }
  }
  sminres_finalize();
  /* 実行結果の検証 */
  double complex *temp;
  temp = (double complex *)calloc(N, sizeof(double complex));
  int ONE = 1;
  for(int k=0; k<M; k++){
    SpMV(A_row,A_col,A_ele, &(x[k*N]), temp, N);
    double complex zTMP = sigma[k];
    zaxpy_(&N, &zTMP, &(x[k*N]), &ONE, temp, &ONE);
    zTMP = -1.0;
    zaxpy_(&N, &zTMP, b, &ONE, temp, &ONE);
    fprintf(stdout, "%d %lf %lf %e\n", k,
            creal(sigma[k]), cimag(sigma[k]),
            dznrm2_(&N, temp, &ONE));
  }
    free(temp);
  return 0;
}

/*
gcc sample2.c ../c/sminres_function.c -lm -lgfortran -lblas -llapack
*/
