#include "definition.h"
#include "function_sminres.h"

void sminres_init()
{
  
}
void sminres_update()
{
  
}
void sminres_finalize()
{
  
}

void solve_sminres(const double complex *sigma, const int M,
                   const int *A_row, const int *A_col, const double complex *A_ele,
                   double complex *x, const double complex *rhs, const int N,
                   int itermax, double threshold, int *status)
{
  const int N2=N*2, N3=N*3, ONE=1;
  double dTMP;         // BLAS用倍精度実数変数
  double complex zTMP; // BLAS用倍精度複素数変数
  int j, k;            // ループ用変数
  // delare variables
  double complex *v;
  double alpha, beta[2];
  double complex *T;
  double *G1;
  double complex *G2;
  double complex *p, *f;
  double *h;
  int conv_num, *is_conv;
  double rhs_nrm;

  // allocate memory
  // 2次元配列から1次元配列にしたが最大値が心配 (ex. M*N3 > INT_MAX)
  v  = (double complex *)calloc(N3,   sizeof(double complex));
  T  = (double complex *)calloc(M*4,   sizeof(double complex));
  G1 = (double         *)calloc(M*3,   sizeof(double));
  G2 = (double complex *)calloc(M*3,   sizeof(double complex));
  p  = (double complex *)calloc(M*N3, sizeof(double complex));
  f  = (double complex *)calloc(M,    sizeof(double complex));
  h  = (double         *)calloc(M,    sizeof(double));
  conv_num = 0;
  is_conv = (int *)calloc(M, sizeof(int));

  // Pre-Process
  rhs_nrm = dznrm2_(&N, rhs, &ONE);
  dTMP = 1 / rhs_nrm;
  zcopy(&N, rhs, &ONE, &(v[N]), &ONE);
  zdscal_(&N, &dTMP, &(v[N]), &ONE);
  beta[0] = 0.0;
  for(k=0; k<M; k++){
    f[k] = 1.0;
    h[k] = rhs_nrm;
  }

  // Main-Process
  for(j=0; j<itermax; j++){
    // Lanczos process
    sample_SpMV(A_row,A_col,A_ele, &(v[N]), &(v[N2]), N);
    zTMP = -beta[0];
    zaxpy_(&N, &zTMP, &(v[0]), &ONE, &(v[N2]), &ONE);
    alpha = creal( zdotc_(&N, &(v[N]), &ONE, &(v[N2]), &ONE) );
    zTMP = -alpha[0];
    zaxpy_(&N, &zTMP, &(v[N]), &ONE, &(v[N2]), &ONE);
    beta[1] = dznrm2_(&N, &(v[N2]), &ONE);
    dTMP = 1 / beta[1];
    zdscal_(&N, &dTMP, &(v[N2]), &ONE);
    // shift systems
    for(k=0; k<M; k++){
      if(is_conv[k] == 1){
        continue;
      }
      T[k*4+0] = 0;
      T[k*4+1] = beta[0]; T[k*4+2] = alpha + sigma[k]; T[k*4+3] = beta[1];
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

      zTMP = rhs_nrm * G1[k*3+2] * f[k];
      zaxpy_(&N, &zTMP, &(p[k*N3+N2]), &ONE, &(x[k*N]), &ONE);

      f[k] = -conj(G2[k*3+2]) * f[k];
      h[k] = cabs(-conj(G2[k*3+2])) * h[k];
      G1[k*3+0]=G1[k*3+1]; G1[k*3+1]=G1[k*3+2];
      G2[k*3+0]=G2[k*3+1]; G2[k*3+1]=G2[k*3+2];

      if(h[k] < threshold){
        conv_num++;
        is_conv[k] = 1;
        continue;
      }
    }
    if(conv_num == M){
      break;
    }
    zcopy_(&N, &(v[N]), &ONE, &(v[0]), &ONE);
    zcopy_(&N, &(v[2*N]), &ONE, &(v[N]), &ONE);
    beta[0] = beta[1];
  }

  // Post-Process
  free(v);
  free(T);
  free(G1); free(G2);
  fre(p); free(f); free(h);
  free(is_conv);
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
