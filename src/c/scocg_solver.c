#include "definition.h"
#include "scocg_solver.h"

int solve_scocg(const double complex *sigma, const int M,
		const int *A_row, const int *A_col, const double complex *A_ele,
		double complex *x, double complex *rhs, int N,
		int itermax, double threshold, int *status)
{
  int ONE=1;
  double complex zTMP; // BLAS用倍精度複素数変数(並列化は考えな
  int j, k;
  // declare variables
  int ret;
  int s, t;
  double complex *r, *p, *Ap, *pi;
  double complex *alpha, *alpha_old, *beta;
  double complex rr, rr_old;
  double *res;
  int conv_num, *is_conv;
  
  // allocate memory
  r         = (double complex *)calloc(M*N, sizeof(double complex));
  p         = (double complex *)calloc(M*N, sizeof(double complex));
  Ap        = (double complex *)calloc(N,   sizeof(double complex));
  pi        = (double complex *)calloc(M*3, sizeof(double complex));
  alpha     = (double complex *)calloc(M,   sizeof(double complex));
  alpha_old = (double complex *)calloc(M,   sizeof(double complex));
  beta      = (double complex *)calloc(M,   sizeof(double complex));
  res       = (double         *)calloc(M,   sizeof(double));
  conv_num  = 0;
  is_conv   = (int            *)calloc(M,   sizeof(int));

  /* shifted COCG method */
  // 変数の初期化
  ret = 0;
  s = 0;
  for(k=0; k<M; k++){
    zcopy_(&N, rhs, &ONE, &(r[k*N]), &ONE);
    pi[k*3+0] = 1;
    pi[k*3+1] = 1;
    alpha[k]  = 1;
    beta[k]   = 0;
  }
  // メインループ
  for(j=0; j<itermax; j++){
    // seed system
    zTMP = beta[s];
    zscal_(&N, &zTMP, &(p[s*N]), &ONE);
    zTMP = 1.0;
    zaxpy_(&N, &zTMP, &(r[s*N]), &ONE, &(p[s*N]), &ONE);

    alpha_old[s] = alpha[s];
    rr = zdotu_(&N, &(r[s*N]), &ONE, &(r[s*N]), &ONE);
    sample_SpMV(A_row,A_col,A_ele, &(p[s*N]), Ap, N);
    zTMP = sigma[s];
    zaxpy_(&N, &zTMP, &(p[s*N]), &ONE, Ap, &ONE);
    alpha[s] = rr / zdotu_(&N, &(p[s*N]), &ONE, Ap, &ONE);

    zTMP = alpha[s];
    zaxpy_(&N, &zTMP, &(p[s*N]), &ONE, &(x[s*N]), &ONE);

    zTMP = -alpha[s];
    zaxpy_(&N, &zTMP, Ap, &ONE, &(r[s*N]), &ONE);
    res[s] = dznrm2_(&N, &(r[s*N]), &ONE);
    if(res[s] < threshold && is_conv[s] != 0){
      conv_num++;
      is_conv[s] = j;
    }
    // add system
    for(k=0; k<M; k++){
      if(is_conv[k] != 0 || k == s){
	continue;
      }
      pi[k*3+2] = (1 + (beta[s]/alpha_old[s])*alpha[s] + alpha[s]*(sigma[k]-sigma[s]))*pi[k*3+1] - (beta[s]/alpha_old[s])*alpha[s]*pi[k*3+0];
      alpha[k] = (pi[k*3+1]/pi[k*3+2])*alpha[s];
      beta[k]  = (pi[k*3+0]*pi[k*3+0])/(pi[k*3+1]*pi[k*3+1])*beta[s];

      zTMP = beta[k];
      zscal_(&N, &zTMP, &(p[k*N]), &ONE);
      zTMP = 1.0;
      zaxpy_(&N, &zTMP, &(r[k*N]), &ONE, &(p[k*N]), &ONE);
      zTMP = alpha[k];
      zaxpy_(&N, &zTMP, &(p[k*N]), &ONE, &(x[k*N]), &ONE);

      zcopy_(&N, &(r[s*N]), &ONE, &(r[k*N]), &ONE);
      zTMP = 1.0 / pi[k*3+2];
      zscal_(&N, &zTMP, &(r[k*N]), &ONE);
      res[k] = dznrm2_(&N, &(r[k*N]), &ONE);

      pi[k*3+0] = pi[k*3+1];
      pi[k*3+1] = pi[k*3+2];

      if(res[k] < threshold){
	conv_num++;
	is_conv[k] = j;
	continue;
      }
    }
    if(conv_num == M){
      ret = 1;
      break;
    }
    rr_old = rr;
    rr = zdotu_(&N, &(r[s*N]), &ONE, &(r[s*N]), &ONE);
    beta[s] = rr / rr_old;
    // seed switching
    if(is_conv[s] == 1){
      t = s;
      for(k=0; k<M; k++) if(res[k] > res[t]) t = k;
      beta[t] = (pi[t*3+0]*pi[t*3+0])/(pi[t*3+1]*pi[t*3+1])*beta[s];
      rr = rr / (pi[t*3+1]*pi[t*3+1]);
      for(k=0; k<M; k++){
	if(k == t) continue;
	pi[k*3+0] = pi[k*3*0] / pi[t*3*0];
	pi[k*3+1] = pi[k*3*1] / pi[t*3*1];
      }
      s = t;
    }
    if(status[1] == 1){
      if(j % status[2] == 0){
	fprintf(stderr, "%d ", j);
        for(k=0; k<M; k++) fprintf(stderr, "%e ", res[k]);
        fprintf(stderr, "\n");
      }
    }
  }
  // 終端処理
  if(status[0] == 1){
    double complex *temp;
    temp = (double complex *)calloc(N, sizeof(double complex));
    if(ret == 1)
      fprintf(stdout, "Converged all equation\n");
    else
      fprintf(stdout, "Unconverged all equation\n");
    for(k=0; k<M; k++){
      sample_SpMV(A_row,A_col,A_ele, &(x[k*N]), temp, N);
      zTMP = sigma[k];
      zaxpy_(&N, &zTMP, &(x[k*N]), &ONE, temp, &ONE);
      zTMP = -1.0;
      zaxpy_(&N, &zTMP, rhs, &ONE, temp, &ONE);
      fprintf(stdout, "%d %d %lf %lf %e %e\n", k, is_conv[k],
              creal(sigma[k]), cimag(sigma[k]),
              res[k], dznrm2_(&N, temp, &ONE));
    }
    free(temp);
  }
  free(r); free(p); free(Ap); free(pi);
  free(alpha); free(alpha_old); free(beta);
  free(res); free(is_conv);
  return ret;
}

void sample_SpMV(const int *A_row, const int *A_col, const double complex *A_ele,
                 const double complex *x, double complex *b, int N)
{
  double complex tmp;
  for(int i=0; i<N; i++){
    tmp = CMPLX(0.0, 0.0);
    for(int j=A_row[i]; j<A_row[i+1]; j++)
      tmp += A_ele[j] * x[A_col[j]];
    b[i] = tmp;
  }
}
