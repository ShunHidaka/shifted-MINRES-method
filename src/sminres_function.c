#include "definition.h"
#include "sminres_function.h"

/* BLAS用定数(並列化は考えない) */
int ONE=1;
double dTMP;
double complex zTMP;

/* 使用する変数 */
int N, M;
double threshold;
double complex *sigma;
double complex *v;
double *alpha, *beta;
double complex *T;
double *G1; double complex *G2;
double complex *p, *f;
double *h, *r0nrm;
int *conv_num, *is_conv;

void sminres_init(const int input_N, double complex *input_rhs,
                  const int input_M, double complex *input_sigma,
                  double complex *input_v, double complex *x,
                  const double input_threshold)
{
  /* 定数の代入 */
  N         = input_N;
  M         = input_M;
  threshold = input_threshold;
  /* メモリの割り当て */
  v     = (double complex *)calloc(N*3,   sizeof(double complex));
  alpha = (double         *)calloc(1,     sizeof(double));
  beta  = (double         *)calloc(2,     sizeof(double));
  T     = (double complex *)calloc(M*4,   sizeof(double complex));
  G1    = (double         *)calloc(M*3,   sizeof(double));
  G2    = (double complex *)calloc(M*3,   sizeof(double complex));
  p     = (double complex *)calloc(M*N*3, sizeof(double complex));
  f     = (double complex *)calloc(M,     sizeof(double complex));
  h     = (double         *)calloc(M,     sizeof(double));
  r0nrm = (double         *)calloc(1,     sizeof(double));
  sigma = (double complex *)calloc(M,     sizeof(double complex));
  conv_num = (int *)calloc(1, sizeof(int));
  is_conv  = (int *)calloc(M, sizeof(int));
  /* 変数の初期化 */
  for(int i=0; i<M*N; i++)
    x[i] = CMPLX(0.0, 0.0);
  r0nrm[0] = dznrm2_(&N, input_rhs, &ONE);
  dTMP = 1 / r0nrm[0];
  zcopy_(&N, input_rhs, &ONE, &(v[N]), &ONE);
  zdscal_(&N, &dTMP, &(v[N]), &ONE);
  beta[0] = 0.0;
  for(int k=0; k<M; k++){
    f[k] = 1.0;
    h[k] = r0nrm[0];
  }
  zcopy_(&M, input_sigma, &ONE, sigma, &ONE);
  /* ユーザー用変数に格納 */
  zcopy_(&N, &(v[N]), &ONE, input_v, &ONE);
}
void sminres_update(int j, double complex *input_v, double complex *input_Av,
                    double complex *x, int *status)
{
  int N2=N*2, N3=N*3;
  // Lanczos過程
  zcopy_(&N, input_Av, &ONE, &(v[2*N]), &ONE);
  zTMP = -beta[0];
  zaxpy_(&N, &zTMP, &(v[0]), &ONE, &(v[N2]), &ONE);
  alpha[0] = creal( zdotc_(&N, &(v[N]), &ONE, &(v[N2]), &ONE) );
  zTMP = -alpha[0];
  zaxpy_(&N, &zTMP, &(v[N]), &ONE, &(v[N2]), &ONE);
  beta[1] = dznrm2_(&N, &(v[N2]), &ONE);
  dTMP = 1 / beta[1];
  zdscal_(&N, &dTMP, &(v[N2]), &ONE);
  // 近似解の更新
  for(int k=0; k<M; k++){
    if(is_conv[k] != 0){
      continue;
    }
    T[k*4+0] = 0;
    T[k*4+1] = beta[0]; T[k*4+2] = alpha[0] + sigma[k]; T[k*4+3] = beta[1];
    if(j >= 2){ zrot_(&ONE, &(T[k*4+0]), &ONE, &(T[k*4+1]), &ONE, &(G1[k*3+0]), &(G2[k*3+0]));}
    if(j >= 1){ zrot_(&ONE, &(T[k*4+1]), &ONE, &(T[k*4+2]), &ONE, &(G1[k*3+1]), &(G2[k*3+1]));}
    zrotg_( &(T[k*4+2]), &(T[k*4+3]), &(G1[k*3+2]), &(G2[k*3+2]) );

    zcopy_(&N, &(p[k*N3+N]) , &ONE, &(p[k*N3+0]) , &ONE);
    zcopy_(&N, &(p[k*N3+N2]), &ONE, &(p[k*N3+N]) , &ONE);
    zcopy_(&N, &(v[N])      , &ONE, &(p[k*N3+N2]), &ONE);
    zTMP = -T[k*4+0];
    zaxpy_(&N, &zTMP, &(p[k*N3+0]), &ONE, &(p[k*N3+N2]), &ONE);
    zTMP = -T[k*4+1];
    zaxpy_(&N, &zTMP, &(p[k*N3+N]), &ONE, &(p[k*N3+N2]), &ONE);
    zTMP = 1 / T[k*4+2];
    zscal_(&N, &zTMP, &(p[k*N3+N2]), &ONE);

    zTMP = r0nrm[0] * G1[k*3+2] * f[k];
    zaxpy_(&N, &zTMP, &(p[k*N3+N2]), &ONE, &(x[k*N]), &ONE);

    f[k] = -conj(G2[k*3+2]) * f[k];
    h[k] =  cabs(-conj(G2[k*3+2])) * h[k];
    G1[k*3+0]=G1[k*3+1]; G1[k*3+1]=G1[k*3+2];
    G2[k*3+0]=G2[k*3+1]; G2[k*3+1]=G2[k*3+2];

    if(h[k] < threshold){
      conv_num[0]++;
      is_conv[k] = 1;
      continue;
    }
  }
  zcopy_(&N, &(v[N]),  &ONE, &(v[0]), &ONE);
  zcopy_(&N, &(v[N2]), &ONE, &(v[N]), &ONE);
  zcopy_(&N, &(v[N2]), &ONE, input_v, &ONE);
  beta[0] = beta[1];
  // 収束判定
  if(conv_num[0] == M){
    *status = 1;
  }
  else{
    *status = 0;
  }
}
void sminres_finalize()
{
  free(v);
  free(alpha); free(beta);
  free(T); free(G1); free(G2);
  free(p); free(f); free(h);
  free(r0nrm);
  free(sigma);
  free(conv_num); free(is_conv);
}
void sminres_getresidual(double *input_h)
{
  zcopy_(&M, &h, &ONE, input_h, &ONE);
}
