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
  int M;
  double complex *sigma;
  set_shifts(&M, &sigma);

  // declare variables and allocate memory
  int s, t;
  double complex **x, **r, **p, *Ap, **pi;
  double complex *alpha, *alpha_old, *beta;
  double complex rr, rr_old;
  double *res;
  // allocate memory
  x         = (double complex **)calloc(M, sizeof(double complex *));
  r         = (double complex **)calloc(M, sizeof(double complex *));
  p         = (double complex **)calloc(M, sizeof(double complex *));
  Ap        = (double complex  *)calloc(N, sizeof(double complex));
  pi        = (double complex **)calloc(M, sizeof(double complex *));
  alpha     = (double complex  *)calloc(M, sizeof(double complex));
  alpha_old = (double complex  *)calloc(M, sizeof(double complex));
  beta      = (double complex  *)calloc(M, sizeof(double complex));
  res       = (double          *)calloc(M, sizeof(double));
  for(k=0; k<M; k++){
    x[k]  = (double complex *)calloc(N, sizeof(double complex));
    r[k]  = (double complex *)calloc(N, sizeof(double complex));
    p[k]  = (double complex *)calloc(N, sizeof(double complex));
    pi[k] = (double complex *)calloc(3, sizeof(double complex));
  }
  int  conv_num = 0;
  int *is_conv  = (int *)calloc(M, sizeof(int));
  double complex *temp = (double complex *)calloc(N, sizeof(double complex));

  // shifted COCG method
  fprintf(stdout, "Please enter seed index (default 0): ");
  if(scanf("%d", &s) != 1 || s >= M){
    s = 0;
  }
  for(k=0; k<M; k++){
    zcopy_(&N, &(b[0]), &ONE, &(r[k][0]), &ONE);
    pi[k][0] = pi[k][1] = 1;
  }
  alpha[s] = 1; beta[s] = 0;
  // Main loop
  time_t start_time, end_time;
  start_time = time(NULL);
  for(j=1; j<MAX_ITR; j++){
    // seed system
    cTMP = beta[s];
    zscal_(&N, &cTMP, &(p[s][0]), &ONE);
    cTMP = 1.0 + 0.0I;
    zaxpy_(&N, &cTMP, &(r[s][0]), &ONE, &(p[s][0]), &ONE);

    alpha_old[s] = alpha[s];
    rr = zdotu_(&N, &(r[s][0]), &ONE, &(r[s][0]), &ONE);
    SpMV(A_row,A_col,A_ele, p[s], Ap, N);
    cTMP = sigma[s];
    zaxpy_(&N, &cTMP, &(p[s][0]), &ONE, &(Ap[0]), &ONE);
    alpha[s] = rr / zdotu_(&N, &(p[s][0]), &ONE, &(Ap[0]), &ONE);

    cTMP = alpha[s];
    zaxpy_(&N, &cTMP, &(p[s][0]),&ONE, &(x[s][0]), &ONE);

    cTMP = -alpha[s];
    zaxpy_(&N, &cTMP, &(Ap[0]), &ONE, &(r[s][0]), &ONE);
    res[s] = dznrm2_(&N, &(r[s][0]), &ONE);

    // shift systems
    for(k=0; k<M; k++){
      if(is_conv[k] != 0 || k == s){
        continue;
      }
      pi[k][2] = (1 + (beta[s]/alpha_old[s])*alpha[s] + alpha[s]*(sigma[k]-sigma[s]))*pi[k][1] - (beta[s]/alpha_old[s])*alpha[s]*pi[k][0];
      alpha[k] = (pi[k][1]/pi[k][2])*alpha[s];
      beta[k] = (pi[k][0]/pi[k][1])*(pi[k][0]/pi[k][1])*beta[s];

      cTMP = beta[k];
      zscal_(&N, &cTMP, &(p[k][0]), &ONE);
      cTMP = 1.0 + 0.0I;
      zaxpy_(&N, &cTMP, &(r[k][0]), &ONE, &(p[k][0]), &ONE);
      cTMP = alpha[k];
      zaxpy_(&N, &cTMP, &(p[k][0]), &ONE, &(x[k][0]), &ONE);

      zcopy_(&N, &(r[s][0]), &ONE, &(r[k][0]), &ONE);
      cTMP = 1.0 / pi[k][2];
      zscal_(&N, &cTMP, &(r[k][0]), &ONE);

      res[k] = dznrm2_(&N, &(r[k][0]), &ONE);
      if(res[k]/r0norm < 1e-13){
        conv_num++;
        is_conv[k]  = j;
      }

      pi[k][0] = pi[k][1];
      pi[k][1] = pi[k][2];
    }
    if(res[s]/r0norm < 1e-13){
      conv_num++;
      is_conv[s]  = j;
    }
    /*
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
    */
    // Determine Convergence
    if(conv_num == M){
      break;
    }
    rr_old = rr;
    rr = zdotu_(&N, &(r[s][0]), &ONE, &(r[s][0]), &ONE);
    beta[s] = rr / rr_old;

    // seed switching
    if(is_conv[s] != 0){
      t = s;
      for(k=0; k<M; k++){
        if(res[k] > res[t] && k != s) t = k;
      }
      beta[t]  = (pi[t][0]/pi[t][1])*(pi[t][0]/pi[t][1])*beta[s];
      for(k=0; k<M; k++){
        if(k == t) continue;
        pi[k][0] = pi[k][0] / pi[t][0];
        pi[k][1] = pi[k][1] / pi[t][1];
      }
      fprintf(stdout, "# SWITCH [%d] to [%d] in %d : %e %e\n", s, t, j, res[s], res[t]);
      s = t;
    }

  }
  end_time = time(NULL);

  // Output results
  fprintf(stdout, "# shifted cocg method with seed switching\n");
  fprintf(stdout, "# matrix=%s\n", FNAME);
  fprintf(stdout, "# iteration=%d, status=%d/%d, time=%ld\n", j, conv_num,M, end_time-start_time);
  fprintf(stdout, "# k, real(sigma[k]), imag(sigma[k]) conv_itr REAL_RES TRUE_RES\n");
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
    free(x[k]); free(r[k]); free(p[k]); free(pi[k]);
  }
  free(x); free(r); free(p); free(Ap);
  free(pi); free(alpha); free(alpha_old);
  free(beta); free(res);
  free(b); free(sigma);
  free(A_row); free(A_col); free(A_ele);
  return 0;
}
