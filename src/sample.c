#include "sminres_functions.h"
#include "function_util.h"

char *FNAME = "VCNT900h_A.csr";

int main(void){
  int m;
  int p1 = 1;
  double complex m1 = -1.0;
  /* read R-Symmetric or C-Hermitian Matrix 'A' */
  int N, DATASIZE;
  int *A_row, *A_col;
  double complex *A_ele;
  read_csr(FNAME, &N, &DATASIZE, &A_row, &A_col, &A_ele);
  /* generate shift points */
  int M=1001;
  double complex *sigma;
  sigma = (double complex *)calloc(M, sizeof(double complex));
  for(m=0; m<M; m++){
    //sigma[m] = 0.01 * cexp( 2 * M_PI * I * (m + 0.5) / M );
    sigma[m] = (-0.501 + 0.001*m) + 0.001I;
  }
  /* set rhs-vector */
  double complex *b;
  b = (double complex *)calloc(N, sizeof(double complex));
  for(int i=0; i<N; i++){
    b[i] = 1.0;
  }
  double bnrm = dznrm2_(&N,b,&p1);

  int itermax = 100000;
  double threshold = 1e-13;
  double complex **x, *q, *Aq;
  double *res;
  int status;
  x   = (double complex **)calloc(M, sizeof(double complex *));
  for(m=0; m<M; m++){
    x[m] = (double complex *)calloc(N, sizeof(double complex));
  }
  q   = (double complex  *)calloc(N, sizeof(double complex));
  Aq  = (double complex  *)calloc(N, sizeof(double complex));
  res = (double *)calloc(M, sizeof(double));
  
  sminres_initialize(N, b, M, sigma, q, threshold);
  for(int j=1; j<itermax; j++){
    SpMV(A_row,A_col,A_ele, q, Aq, N);
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
  for(m=0; m<M; m++){
    SpMV(A_row,A_col,A_ele, &(x[m][0]), TMP, N);
    zaxpy_(&N, &(sigma[m]), &(x[m][0]), &p1, TMP, &p1);
    zaxpy_(&N, &m1, b, &p1, TMP, &p1);
    fprintf(stdout, "%d %lf %lf %e %e\n",
            m, creal(sigma[m]), cimag(sigma[m]), res[m]/bnrm, dznrm2_(&N,TMP,&p1)/bnrm);
  }
  return 0;
}
