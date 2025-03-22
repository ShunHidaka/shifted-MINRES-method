#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include "function_util.h"
#include "function_blas.h"

int main(int argc, char *argv[])
{
  int j, k;
  // prepare Matrix A
  int N, DATASIZE;
  int *A_row, *A_col; double complex *A_ele;
  set_fname(argc, argv);
  read_csr(FNAME, &N, &DATASIZE, &A_row, &A_col, &A_ele);
  // prepare right-hand side vector b
  double complex *b;
  double r0norm;
  set_rhsVector(N, &b, &r0norm);
  // prepare shifts sigma
  int M=1;
  double complex *sigma;
  set_shifts(&M, &sigma);

  // delcare variables and allocate memory
  int s=0;
  double complex **x, **r1, **p1, **pi;
  double complex *alpha, *alpha_old, *beta;
  double complex *r2, **p2, *Ap1, *Ap2;
  double complex r1r2, r1r2_old;
  double *res;
  x         = (double complex **)calloc(M, sizeof(double complex *));
  r1        = (double complex **)calloc(M, sizeof(double complex *));
  r2        = (double complex  *)calloc(N, sizeof(double complex));
  p1        = (double complex **)calloc(M, sizeof(double complex *));
  p2        = (double complex **)calloc(M, sizeof(double complex *));
  Ap1       = (double complex  *)calloc(N, sizeof(double complex));
  Ap2       = (double complex  *)calloc(N, sizeof(double complex));
  pi        = (double complex **)calloc(M, sizeof(double complex));
  alpha     = (double complex  *)calloc(M, sizeof(double complex));
  alpha_old = (double complex  *)calloc(M, sizeof(double complex));
  beta      = (double complex  *)calloc(M, sizeof(double complex));
  res       = (double          *)calloc(M, sizeof(double));
  for(k=0; k<M; k++){
    x[k]  = (double complex *)calloc(N, sizeof(double complex));
    r1[k] = (double complex *)calloc(N, sizeof(double complex));
    p1[k] = (double complex *)calloc(N, sizeof(double complex));
    p2[k] = (double complex *)calloc(N, sizeof(double complex));
    pi[k] = (double complex *)calloc(3, sizeof(double complex));
  }
  int conv_num = 0;
  int *is_conv = (int *)calloc(M, sizeof(int));
  double complex *temp = (double complex *)calloc(N, sizeof(double complex));

  // shifted BiCG method
  s = 4;
  for(k=0; k<M; k++){
    zcopy_(&N, &(b[0]), &ONE, &(r1[k][0]), &ONE);
    pi[k][0] = pi[k][1] = 1;
  }
  zcopy_(&N, &(r1[s][0]), &ONE, &(r2[0]), &ONE);
  alpha[s] = 1; beta[s] = 0;
  // Main loop
  time_t start_time, end_time;
  start_time = time(NULL);
  for(j=1; j<MAX_ITR; j++){
    // seed system
    cTMP = beta[s];
    zscal_(&N, &cTMP, &(p1[s][0]), &ONE);
    cTMP = 1.0 + 0.0I;
    zaxpy_(&N, &cTMP, &(r1[s][0]), &ONE, &(p1[s][0]), &ONE);
    cTMP = conj(beta[s]);
    zscal_(&N, &cTMP, &(p2[s][0]), &ONE);
    cTMP = 1.0 + 0.0I;
    zaxpy_(&N, &cTMP, &(r2[0]), &ONE, &(p2[s][0]), &ONE);

    alpha_old[s] = alpha[s];
    r1r2 = zdotc_(&N, &(r2[0]), &ONE, &(r1[s][0]), &ONE); // <r1[s], r2> = r2^H * r1[s]
    SpMV(A_row,A_col,A_ele, p1[s], Ap1, N);
    cTMP = sigma[s];
    zaxpy_(&N, &cTMP, &(p1[s][0]), &ONE, &(Ap1[0]), &ONE);
    alpha[s] = r1r2 / zdotc_(&N, &(p2[s][0]), &ONE, &(Ap1[0]), &ONE);

    cTMP = alpha[s];
    zaxpy_(&N, &cTMP, &(p1[s][0]), &ONE, &(x[s][0]), &ONE);

    cTMP = -alpha[s];
    zaxpy_(&N, &cTMP, &(Ap1[0]), &ONE, &(r1[s][0]), &ONE);
    res[s] = dznrm2_(&N, &(r1[s][0]), &ONE);

    // shift systems
    for(k=0; k<M; k++){
      if(is_conv[k] != 0 || k == s){
        continue;
      }
      pi[k][2] = (1 + (beta[s]/alpha_old[s])*alpha[s] + alpha[s]*(sigma[k]-sigma[s]))*pi[k][1] - (beta[s]/alpha_old[s])*alpha[s]*pi[k][0];
      alpha[k] = (pi[k][1]/pi[k][2])*alpha[s];
      beta[k] = (pi[k][0]/pi[k][1])*(pi[k][0]/pi[k][1])*beta[s];
      // p1[k]
      cTMP = beta[k];
      zscal_(&N, &cTMP, &(p1[k][0]), &ONE);
      cTMP = 1.0 + 0.0I;
      zaxpy_(&N, &cTMP, &(r1[k][0]), &ONE, &(p1[k][0]), &ONE);
      // p2[k]
      cTMP = conj(beta[k]);
      zscal_(&N, &cTMP, &(p2[k][0]), &ONE);
      cTMP = 1.0 / conj(pi[k][1]);
      zaxpy_(&N, &cTMP, &(r2[0]), &ONE, &(p2[k][0]), &ONE);
      // x[k]
      cTMP = alpha[k];
      zaxpy_(&N, &cTMP, &(p1[k][0]), &ONE, &(x[k][0]), &ONE);

      zcopy_(&N, &(r1[s][0]), &ONE, &(r1[k][0]), &ONE);
      cTMP = 1.0 / pi[k][2];
      zscal_(&N, &cTMP, &(r1[k][0]), &ONE);

      res[k] = dznrm2_(&N, &(r1[k][0]), &ONE);
      if(res[k]/r0norm < 1e-13){
        conv_num++;
        is_conv[k] = j;
      }

      pi[k][0] = pi[k][1];
      pi[k][1] = pi[k][2];
    }
    if(res[s]/r0norm < 1e-13 && is_conv[s] == 0){
      conv_num++;
      is_conv[s]  = j;
    }

    if(j % 5 == 1){
      fprintf(stderr, "%d", j);
      for(k=0; k<M; k++){
        SpMV(A_row,A_col,A_ele, &(x[k][0]), temp, N);
        cTMP = sigma[k];
        zaxpy_(&N, &cTMP, &(x[k][0]), &ONE, temp, &ONE);
        cTMP = -1.0;
        zaxpy_(&N, &cTMP, b, &ONE, temp, &ONE);
        fprintf(stderr, " %e %e", res[k]/r0norm, dznrm2_(&N,temp,&ONE)/r0norm);
      }
      fprintf(stderr, "\n");
    }

    // Determin Convergence
    if(conv_num == M){
      break;
    }
    SpMV(A_row,A_col,A_ele, p2[s], Ap2, N);
    cTMP = conj(sigma[s]);
    zaxpy_(&N, &cTMP, &(p2[s][0]), &ONE, &(Ap2[0]), &ONE);
    cTMP = -conj(alpha[s]);
    zaxpy_(&N, &cTMP, &(Ap2[0]), &ONE, &(r2[0]), &ONE);
    r1r2_old = r1r2;
    r1r2 = zdotc_(&N, &(r2[0]), &ONE, &(r1[s][0]), &ONE);
    beta[s] = r1r2 / r1r2_old;
  }
  end_time = time(NULL);
  
  // Output results
  fprintf(stdout, "# shifted bicg method (seed=%d)\n", s);
  fprintf(stdout, "# matrix=%s\n", FNAME);
  fprintf(stdout, "# iteration=%d, status=%d/%d, time=%ld\n", j, conv_num, M, end_time-start_time);
  for(k=0; k<M; k++){
    SpMV(A_row,A_col,A_ele, &(x[k][0]), temp, N);
    cTMP = sigma[k];
    zaxpy_(&N, &cTMP, &(x[k][0]), &ONE, temp, &ONE);
    cTMP = -1.0;
    zaxpy_(&N, &cTMP, b, &ONE, temp, &ONE);
    fprintf(stdout, "%d %lf %lf %d %e %e\n",
            k, creal(sigma[k]),cimag(sigma[k]), is_conv[k], res[k]/r0norm, dznrm2_(&N,temp,&ONE)/r0norm);
  }
  free(temp);

  free(is_conv);
  for(k=0; k<M; k++){
    free(x[k]); free(r1[k]);
    free(p1[k]); free(p2[k]); free(pi[k]);
  }
  free(x); free(r1); free(r2);
  free(p1); free(p2); free(pi);
  free(alpha); free(alpha_old);
  free(beta); free(Ap1); free(Ap2);
  free(res); free(sigma); free(b);
  free(A_row); free(A_col); free(A_ele);
  return 0;
}
