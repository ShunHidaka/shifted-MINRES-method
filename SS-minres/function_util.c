#include "definition.h"
#include "function_util.h"

void solve_sminres(double complex *sigma,
                   int *A_row, int *A_col, double complex *A_ele,
                   double complex *rhs,
                   double complex *x, int M, int N,
                   int itermax, double threshold, int *status)
{
  int j, k;
  double complex *v;
  double *alpha, *beta;
  double complex **T;
  double **G1;
  double complex **G2;
  double complex  **p, *f;
  double *h, *res;
  int conv_num, *is_conv;
  // allocate memory
  v     = (double complex *)calloc(N*3, sizeof(double complex));
  alpha = (double *)calloc(1, sizeof(double));
  beta  = (double *)calloc(2, sizeof(double));
  T     = (double complex **)calloc(M, sizeof(double complex *));
  G1    = (double **)calloc(M, sizeof(double *));
  G2    = (double complex **)calloc(M, sizeof(double complex *));
  p     = (double complex **)calloc(M, sizeof(double complex *));
  f     = (double complex *)calloc(M, sizeof(double complex));
  h     = (double *)calloc(M, sizeof(double));
  res   = (double *)calloc(M, sizeof(double));
  for(k=0; k<M; k++){
    T[k]  = (double complex *)calloc(4, sizeof(double complex));
    G1[k] = (double *)calloc(3, sizeof(double));
    G2[k] = (double complex *)calloc(3, sizeof(double complex));
    p[k]  = (double complex *)calloc(N*3, sizeof(double complex));
  }
  is_conv = (int *)calloc(M, sizeof(int));

  int ONE = 1;
  double dTMP;
  double complex cTMP;
  // Pre-Process
  double rhs_nrm = dznrm2_(&N, rhs, &ONE);
  dTMP = 1 / rhs_nrm;
  zcopy_(&N, &(rhs[0]), &ONE, &(v[N]), &ONE);
  zdscal_(&N, &dTMP, &(v[N]), &ONE);
  beta[0] = 0.0;
  conv_num = 0;
  for(k=0; k<M; k++){
    f[k] = 1;
    h[k] = rhs_nrm;
    is_conv[k] = 0;
  }
  // Main-Process
  for(j=0; j<itermax; j++){
    // Lanczos process
    SpMV(A_row,A_col,A_ele, &(v[N]), &(v[2*N]), N);
    cTMP = -beta[0];
    zaxpy_(&N, &cTMP, &(v[0]), &ONE, &(v[2*N]), &ONE);
    alpha[0] = creal( zdotc_(&N, &(v[N]), &ONE, &(v[2*N]), &ONE) );
    cTMP = -alpha[0];
    zaxpy_(&N, &cTMP, &(v[N]), &ONE, &(v[2*N]), &ONE);
    beta[1] = dznrm2_(&N, &(v[2*N]), &ONE);
    dTMP = 1 / beta[1];
    zdscal_(&N, &dTMP, &(v[2*N]), &ONE);
    // shift systems
    for(k=0; k<M; k++){
      if(is_conv[k] == 1) continue;

      T[k][0] = 0;
      T[k][1] = beta[0]; T[k][2] = alpha[0] + sigma[k]; T[k][3] = beta[1];
      if(j >= 2){ zrot_(&ONE, &(T[k][0]), &ONE, &(T[k][1]), &ONE, &(G1[k][0]), &(G2[k][0]));}
      if(j >= 1){ zrot_(&ONE, &(T[k][1]), &ONE, &(T[k][2]), &ONE, &(G1[k][1]), &(G2[k][1]));}
      zrotg_( &(T[k][2]), &(T[k][3]), &(G1[k][2]), &(G2[k][2]) );

      zcopy_(&N, &(p[k][N]),   &ONE, &(p[k][0]),   &ONE);
      zcopy_(&N, &(p[k][2*N]), &ONE, &(p[k][N]),   &ONE);
      zcopy_(&N, &(v[N]),      &ONE, &(p[k][2*N]), &ONE);
      cTMP = -T[k][0];
      zaxpy_(&N, &cTMP, &(p[k][0]), &ONE, &(p[k][2*N]), &ONE);
      cTMP = -T[k][1];
      zaxpy_(&N, &cTMP, &(p[k][N]), &ONE, &(p[k][2*N]), &ONE);

      cTMP = 1 / T[k][2];
      zscal_(&N, &cTMP, &(p[k][2*N]), &ONE);
      cTMP = rhs_nrm * G1[k][2] * f[k];
    //cTMP = rhs_nrm * G1[k][2] * f[k] / T[k][2];
      zaxpy_(&N, &cTMP, &(p[k][2*N]), &ONE, &(x[k*N]), &ONE);
      
      f[k] = -conj(G2[k][2])*f[k];
      h[k]    = cabs( -conj(G2[k][2]) ) * h[k];
      G1[k][0]=G1[k][1]; G1[k][1]=G1[k][2];
      G2[k][0]=G2[k][1]; G2[k][1]=G2[k][2];

      if(h[k] < threshold){
        conv_num++;
        is_conv[k] = 1;
        continue;
      }
    }
    // Determine Convergence
    if(conv_num == M){
      status[0] = 1;
      status[1] = j;
      break;
    }
    zcopy_(&N, &(v[N]), &ONE, &(v[0]), &ONE);
    zcopy_(&N, &(v[2*N]), &ONE, &(v[N]), &ONE);
    beta[0] = beta[1];
  }
  // Post-Process
  //  if(status[2] == 1){

  // deallocate memory
  free(v);
  free(alpha); free(beta);
  for(k=0; k<M; k++){
    free(T[k]);
    free(G1[k]); free(G2[k]);
    free(p[k]);
  }
  free(T);free(G1);free(G2);
  free(p);free(f);free(h);
}

void SpMV(const int *A_row, const int *A_col, const double complex *A_ele,
          const double complex *x, double complex *b, int N)
{
  double complex tmp;
  for(int i=0; i<N; i++){
    tmp = 0.0+0.0I;
    for(int j=A_row[i]; j<A_row[i+1]; j++)
      tmp += A_ele[j]*x[A_col[j]];
    b[i] = tmp;
  }
}
FILE *fopen_mtx(const char *fname, const char *mode, int *row_size, int *col_size, int *ele_size)
{
  FILE *fp = NULL;
  fp = fopen(fname, mode);
  if(fp == NULL){
    fprintf(stderr, "Can not open file : %s\n", fname);
    exit(1);
  }
  char chr;
  while( (chr = fgetc(fp)) != EOF && (chr=='%' || chr=='#') )
    while( (chr = fgetc(fp)) != EOF )
      if(chr == '\n') break;
  fseek(fp, -sizeof(char), SEEK_CUR);
  int num, tmp1, tmp2, tmp3;
  num = fscanf(fp, "%d %d %d", &tmp1, &tmp2, &tmp3);
  if(num != 3){ fprintf(stderr, "fopen err\n"); exit(1);}
  if(col_size != NULL) *row_size = tmp1;
  if(row_size != NULL) *col_size = tmp2;
  if(ele_size != NULL) *ele_size = tmp3;
  return fp;
}
void read_csr(const char *fname, int *N, int *DATASIZE,
              int **row_ptr, int **col_ind, double complex **element)
{
  FILE *fp;
  int rsize, csize, esize;
  fp = fopen_mtx(fname, "r", &rsize, &csize, &esize);
  *N = rsize-1; *DATASIZE = esize;

  *row_ptr = (int *)calloc(rsize, sizeof(int));
  *col_ind = (int *)calloc(csize, sizeof(int));
  *element = (double complex *)calloc(esize, sizeof(double complex));

  int row, col, i=0;
  double real, imag=0;
  while( fscanf(fp, "%d %d %lf %lf", &row, &col, &real, &imag) != EOF ){
    if(i < rsize) (*row_ptr)[i] = row;
    if(i < csize) (*col_ind)[i] = col;
    if(i < esize) (*element)[i] = real + imag*I;
    i = i + 1;
  }
  fclose(fp);
}
