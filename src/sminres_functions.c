#include "definition.h"
#include "sminres_functions.h"

/* BLAS用定数(並列化は考えない) */
int ONE=1;
double dTMP;
double complex cTMP;

/* 使用する変数 */
int j, m, N, M;
int N2, N3;
double threshold;
double complex *sigma;
double complex *q;
double *alpha, *beta;
double complex **r;
double **c;
double complex **s;
double complex **p;
double complex *f;
double *h;
double r0nrm;
int conv_num;
int *is_conv;

void sminres_initialize(const int input_N, double complex *input_rhs,
                        const int input_M, double complex *input_sigma,
                        double complex *input_q,
                        const double input_threshold)
{
  /* 定数の代入 */
  N         = input_N;
  M         = input_M;
  N2        = 2*N;
  N3        = 3*N;
  threshold = input_threshold;
  /* メモリの割り当て */
  sigma = (double complex  *)calloc(M,     sizeof(double complex));
  q     = (double complex  *)calloc(N*3, sizeof(double complex));
  alpha = (double          *)calloc(1,   sizeof(double));
  beta  = (double          *)calloc(2,   sizeof(double));
  r     = (double complex **)calloc(M,   sizeof(double complex));
  c     = (double         **)calloc(M,   sizeof(double));
  s     = (double complex **)calloc(M,   sizeof(double complex));
  p     = (double complex **)calloc(M,   sizeof(double complex));
  f     = (double complex  *)calloc(M,   sizeof(double complex));
  h     = (double          *)calloc(M,   sizeof(double));
  for(m=0; m<M; m++){
    r[m] = (double complex *)calloc(3,   sizeof(double complex));
    c[m] = (double         *)calloc(3,   sizeof(double));
    s[m] = (double complex *)calloc(3,   sizeof(double complex));
    p[m] = (double complex *)calloc(N*3, sizeof(double complex));
  }
  is_conv  = (int *)calloc(M, sizeof(int));
  /* 変数の初期化 */
  conv_num = 0;
  r0nrm = dznrm2_(&N, input_rhs, &ONE);
  dTMP = 1 / r0nrm;
  zcopy_(&N, input_rhs, &ONE, &(q[N]), &ONE);
  zdscal_(&N, &dTMP, &(q[N]), &ONE);
  beta[0] = 0.0;
  for(m=0; m<M; m++){
    f[m] = 1.0;
    h[m] = r0nrm;
  }
  zcopy_(&M, input_sigma, &ONE, sigma, &ONE);
  /* ユーザー用変数に格納 */
  zcopy_(&N, &(q[N]), &ONE, input_q, &ONE);
}
//
// j=1; j<MAX_ITR; j++ なら j>=3, j>=2
//
void sminres_update(int j, double complex *input_q, double complex *input_Aq,
                    double complex **x, int *status)
{
  // Lanczos過程
  zcopy_(&N, input_Aq, &ONE, &(q[2*N]), &ONE);
  cTMP = -beta[0];
  zaxpy_(&N, &cTMP, &(q[0]), &ONE, &(q[N2]), &ONE);
  alpha[0] = creal( zdotc_(&N, &(q[N]), &ONE, &(q[N2]), &ONE) );
  cTMP = -alpha[0];
  zaxpy_(&N, &cTMP, &(q[N]), &ONE, &(q[N2]), &ONE);
  beta[1] = dznrm2_(&N, &(q[N2]), &ONE);
  // 近似解の更新
  for(m=0; m<M; m++){
    if(is_conv[m] != 0){
      continue;
    }
    r[m][0]=0; r[m][1]=beta[0]; r[m][2]=alpha[0]+sigma[m];
    if(j >= 3){ zrot_(&ONE, &(r[m][0]), &ONE, &(r[m][1]), &ONE, &(c[m][0]), &(s[m][0]));}
    if(j >= 2){ zrot_(&ONE, &(r[m][1]), &ONE, &(r[m][2]), &ONE, &(c[m][1]), &(s[m][1]));}
    cTMP = beta[1];
    zrotg_(&(r[m][2]), &cTMP, &(c[m][2]), &(s[m][2]));

    zcopy_(&N, &(p[m][N]),   &ONE, &(p[m][0]),   &ONE);
    zcopy_(&N, &(p[m][N2]), &ONE, &(p[m][N]),   &ONE);
    zcopy_(&N, &(q[N]),      &ONE, &(p[m][N2]), &ONE);
    cTMP = -r[m][0];
    zaxpy_(&N, &cTMP, &(p[m][0]), &ONE, &(p[m][N2]), &ONE);
    cTMP = -r[m][1];
    zaxpy_(&N, &cTMP, &(p[m][N]), &ONE, &(p[m][N2]), &ONE);
    cTMP = 1.0 / r[m][2];
    zscal_(&N, &cTMP, &(p[m][N2]), &ONE);
    cTMP = r0nrm*c[m][2]*f[m];
    zaxpy_(&N, &cTMP, &(p[m][N2]), &ONE, &(x[m][0]), &ONE);

    f[m] = -conj(s[m][2])*f[m];
    h[m] = cabs(-conj(s[m][2]))*h[m];
    c[m][0]=c[m][1]; c[m][1]=c[m][2];
    s[m][0]=s[m][1]; s[m][1]=s[m][2];

    if(h[m]/r0nrm < threshold){
      conv_num++;
      is_conv[m] = 1;
      continue;
    }
  }

  dTMP = 1 / beta[1];
  zdscal_(&N, &dTMP, &(q[N2]), &ONE);
  beta[0] = beta[1];
  zcopy_(&N, &(q[N]),  &ONE, &(q[0]), &ONE);
  zcopy_(&N, &(q[N2]), &ONE, &(q[N]), &ONE);
  zcopy_(&N, &(q[N2]), &ONE, input_q, &ONE);

  // 収束判定
  if(conv_num == M){
    *status = 1;
  }
  else{
    *status = 0;
  }
}
void sminres_finalize()
{
  free(sigma);
  free(q);
  free(alpha); free(beta);
  for(m=0; m<M; m++){
    free(r[m]);
    free(c[m]);
    free(s[m]);
    free(p[m]);
  }
  free(r); free(c); free(s);
  free(p);
  free(f); free(h);
  free(is_conv);
}
void sminres_getresidual(double *input_h)
{
  dcopy_(&M, h, &ONE, input_h, &ONE);
}
