#include "sminres_functions.h"
#include "function_util.h"

//char *FNAME = "ELSES_MATRIX_DIAB18h_A.mtx";
char *FNAME = "ELSES_MATRIX_VCNT400std_A.mtx";
void zhpmv_(char *uplo, int *n, double complex *alpha, double complex *A,
            double complex *x, int *incx, double complex *beta, double complex *y, int *incy);

int main(void){
  /* BLAS CONSTANTS */
  int BLAS_INT_P1 = 1;
  double complex BLAS_CMPLX_P1 =  1.0;
  double complex BLAS_CMPLX_0  =  0.0;
  double complex BLAS_CMPLX_M1 = -1.0;

  /* read R-Symmetric or C-Hermitian Matrix 'A' */
  int N;
  double complex *A = NULL;
  read_mtx(FNAME, &N, &A);
  /* generate shift points */
  int M=10;
  double complex *sigma = NULL;
  set_shifts(&M, &sigma);
  /* set rhs-vector */
  double bnrm;
  double complex *b = NULL;
  set_rhs(N, &b);
  bnrm = dznrm2_(&N, b, &BLAS_INT_P1);

  /* Prepare shifted MINRES method */
  int itermax = 100000;
  double threshold = 1e-13;
  double complex **x, *q, *Aq;
  double *res;
  int status;
  x   = (double complex **)calloc(M, sizeof(double complex *));
  for(int m=0; m<M; m++){
    x[m] = (double complex *)calloc(N, sizeof(double complex));
  }
  q   = (double complex  *)calloc(N, sizeof(double complex));
  Aq  = (double complex  *)calloc(N, sizeof(double complex));
  res = (double *)calloc(M, sizeof(double));

  /* Solve shifted linear systems by shifted MINRES method */
  sminres_initialize(N, b, M, sigma, q, threshold);
  for(int j=1; j<itermax; j++){
    zhpmv_("U", &N, &BLAS_CMPLX_P1, A, q, &BLAS_INT_P1, &BLAS_CMPLX_0, Aq, &BLAS_INT_P1);
    sminres_update(j, q, Aq, x, &status);
    if(status == 1){
      printf("Converged %d\n", j);
      break;
    }
  }
  sminres_getresidual(res);
  sminres_finalize();

  /* verificate approximate solutions */
  double complex *TMP = (double complex *)calloc(N, sizeof(double complex));
  for(int m=0; m<M; m++){
    zhpmv_("U", &N, &BLAS_CMPLX_P1, A, x[m], &BLAS_INT_P1, &BLAS_CMPLX_0, TMP, &BLAS_INT_P1);
    zaxpy_(&N, &(sigma[m]), &(x[m][0]), &BLAS_INT_P1, TMP, &BLAS_INT_P1);
    zaxpy_(&N, &BLAS_CMPLX_M1, b, &BLAS_INT_P1, TMP, &BLAS_INT_P1);
    fprintf(stdout, "%d %lf %lf %e %e\n",
            m, creal(sigma[m]), cimag(sigma[m]), res[m]/bnrm, dznrm2_(&N,TMP,&BLAS_INT_P1)/bnrm);
  }
  return 0;
}
