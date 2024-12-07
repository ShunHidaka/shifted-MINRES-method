#include "definition.h"
#include "sbicg_solver.h"

int solve_sbicg(const double complex *sigma, const int M,
		const int *A_row, const int *A_col, const double complex *A_ele,
		double complex *x, double complex *rhs, int N,
		int itermax, double threshold, int *status)
{
  int ONE=1;
  double complex zTMP; // BLAS用倍精度複素数変数(並列化は考えていない)
  int j, k;
  // declare variables
  int ret;
  int s, t;
  double complex *r1, *r2, *p1, *p2, *Ap1, *Ap2, *pi;
  double complex *alpha, *alpha_old, *beta;
  double complex r1r2, r1r2_old;
  double *res;
  int conv_num, *is_conv;

  // allocate memory
  r1        = (double complex *)calloc(M*N, sizeof(double complex));
  p1        = (double complex *)calloc(M*N, sizeof(double complex));
  Ap1       = (double complex *)calloc(N,   sizeof(double complex));
  r2        = (double complex *)calloc(N,   sizeof(double complex));
  p2        = (double complex *)calloc(M*N, sizeof(double complex));
  Ap2       = (double complex *)calloc(N,   sizeof(double complex));
  alpha     = (double complex *)calloc(M,   sizeof(double complex));
  alpha_old = (double complex *)calloc(M,   sizeof(double complex));
  beta      = (double complex *)calloc(M,   sizeof(double complex));
  res       = (double         *)calloc(M,   sizeof(double));
  is_conv   = (int            *)calloc(M,   sizeof(int));

  /* shifted Bi-CG method*/
  // 変数の初期化
  ret = 0;
  conv_num = 0;
  s = 0;
  for(k=0; k<M; k++){
    zcopy_(&N, rhs, &ONE, &(r1[k*N]), &ONE);
    pi[k*3+0] = 1;
    pi[k*3+1] = 1;
    alpha[k]  = 1;
    beta[k]   = 0;
  }
  zcopy_(&N, &(r1[s*N]), &ONE, &(r2[0]), &ONE);
  // メインループ
  for(j=0; j<itermax; j++){
    // seed system
    zTMP = beta[s];
    zscal_(&N, &zTMP, &(p1[s*N]), &ONE);
    zTMP = 1.0;
    zaxpy_(&N, &zTMP, &(r1[s*N]), &ONE, &(p1[s*N]), &ONE);
    zTMP = conj(beta[s]);
    zscal_(&N, &zTMP, &(p2[s*N]), &ONE);
    zTMP = 1.0;
    zaxpy_(&N, &zTMP, &(r2[s*N]), &ONE, &(p2[s*N]), &ONE);

    alpha_old[s] = alpha[s];
    r1r2 = zdotc_(&N, &(r2[0]), &ONE, &(r1[s][0]), &ONE);
    sample_SpMV(A_row,A_col,A_ele, &(p1[s*N]), Ap1, &ONE);
    zTMP = sigma[s];
    zaxpy_(&N, &cTMP, &(p1[s][0]), &ONE, &(Ap1[0]), &ONE);

    zTMP = alphs[s];
    zaxpy_(&N, &zTMP, &(p1[s*N]), &ONE, &(x[s*N]), &ONE);

    zTMP = -alpha[s];
    zaxpy_(&N, &zTMP, &(Ap1[0]), &ONE, &(x[s*N]), &ONE);
    res[s] = dznrm2_(&N, &(r1[s*N]), &ONE);
    res[s] = log10(res[s] / r0norm);
    if(res[s] < threshold){
      conv_num++;
      is_conv[s] = j;
    }
    // add system
    for(k=0; k<M; k++){
      if(is_conv[k] != 0 || k == s){
	continue;
      }
      
    }
    // seed switching
    if(is_conv[s] != 0){
      t = s;
      for(k=0; k<M; k++) if(res[k] > res[t]) t = k;
      beta[t] = (pi[t*3+0]*pi[t*3+0])/(pi[t*3+1]*pi[t*3+1])*beta[s];
      zTMP = 1.0 / conj(pi[t*3+1]);
      zscal_(&N, &zTMP, &(r2[0]), &ONE);
      r1r2 = r1r2 / (conj(pi[t*3+1])*pi[t*3+1]);
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
