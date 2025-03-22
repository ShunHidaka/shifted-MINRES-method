#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include "function_util.h"
#include "function_blas.h"

int main(int argc, char *argv[])
{
  int j, m;
  // prepare Matrix A
  int N, DATASIZE;
  int *A_row, *A_col; double complex *A_ele;
  set_fname(argc, argv);
  read_csr(FNAME, &N, &DATASIZE, &A_row, &A_col, &A_ele);
  // prepare rhs vector b
  double complex *b;
  double bnrm;
  set_rhsVector(N, &b, &bnrm);
  // prepare shifts \sigma
  int M=1;
  double complex *sigma;
  set_shifts(&M, &sigma);

  // declare variables and allocate memory
  double complex *q;             // q[0]=q_{j-1}, q[N]=q_{j}, q[2*N]=q_{j+1}
  double         alpha, beta[2]; // beta[0]=beta_{j-1}, beta[1]=beta_{j}
  double complex **r;            // r[m][0]=r_{j-2,j}^{(m)}, r[m][1]=r_{j-1,j}^{(m)}, r[m][2]=r_{j,j}^{(m)}
  double         **c;            // element of Givens matrix
  double complex **s;            // {c[m][0, 1, 2],s[m][0, 1, 2]}=G_{j-2}^{(m)}, G_{j-1}^{(m)}, G_{j}^{(m)}
  double complex **p, **x, *f;   // p[m][0], p[m][N], p[m][2*N] = p_{j-2}^{(m)}, p_{j-1}^{(m)}, p_{j}^{(m)}
  double         *h;             //
  q     = (double complex  *)calloc(N*3, sizeof(double complex));
  r     = (double complex **)calloc(M,   sizeof(double complex));
  c     = (double         **)calloc(M,   sizeof(double *));
  s     = (double complex **)calloc(M,   sizeof(double complex *));
  p     = (double complex **)calloc(M,   sizeof(double complex *));
  x     = (double complex **)calloc(M,   sizeof(double complex *));
  f     = (double complex  *)calloc(M,   sizeof(double complex));
  h     = (double          *)calloc(M,   sizeof(double));
  for(m=0; m<M; m++){
    r[m] = (double complex *)calloc(3,   sizeof(double complex));
    c[m] = (double         *)calloc(3,   sizeof(double));
    s[m] = (double complex *)calloc(3,   sizeof(double complex));
    p[m] = (double complex *)calloc(N*3, sizeof(double complex));
    x[m] = (double complex *)calloc(N,   sizeof(double complex));
  }
  int conv_num = 0;
  int *is_conv = (int *)calloc(M, sizeof(int));
  double complex *res = (double complex *)calloc(N, sizeof(double complex));

  // Initialize
  dTMP = 1.0 / bnrm;
  zcopy_(&N, b, &ONE, &(q[N]), &ONE);
  zdscal_(&N, &dTMP, &(q[N]), &ONE);
  beta[0] = 0.0;
  for(m=0; m<M; m++){
    f[m] = 1.0;
    h[m] = bnrm;
  }
  // Main Loop
  time_t start_time, end_time;
  start_time = time(NULL);
  for(j=1; j<MAX_ITR; j++){
    SpMV(A_row,A_col,A_ele, &(q[N]), &(q[2*N]), N);
    cTMP = -beta[0];
    zaxpy_(&N, &cTMP, &(q[0]), &ONE, &(q[2*N]), &ONE);
    alpha = creal( zdotc_(&N, &(q[2*N]), &ONE, &(q[N]), &ONE) );
    cTMP = -alpha;
    zaxpy_(&N, &cTMP, &(q[N]), &ONE, &(q[2*N]), &ONE);
    beta[1] = dznrm2_(&N, &(q[2*N]), &ONE);
    for(m=0; m<M; m++){
      if(is_conv[m] != 0){
        continue;
      }
      r[m][0]=0; r[m][1]=beta[0]; r[m][2]=alpha+sigma[m];
      if(j >= 3){ zrot_(&ONE, &(r[m][0]), &ONE, &(r[m][1]), &ONE, &(c[m][0]), &(s[m][0]));}
      if(j >= 2){ zrot_(&ONE, &(r[m][1]), &ONE, &(r[m][2]), &ONE, &(c[m][1]), &(s[m][1]));}
      cTMP = beta[1];
      zrotg_(&(r[m][2]), &cTMP, &(c[m][2]), &(s[m][2]));
      zcopy_(&N, &(p[m][N]),   &ONE, &(p[m][0]),   &ONE);
      zcopy_(&N, &(p[m][2*N]), &ONE, &(p[m][N]),   &ONE);
      zcopy_(&N, &(q[N]),      &ONE, &(p[m][2*N]), &ONE);
      cTMP = -r[m][0];
      zaxpy_(&N, &cTMP, &(p[m][0]), &ONE, &(p[m][2*N]), &ONE);
      cTMP = -r[m][1];
      zaxpy_(&N, &cTMP, &(p[m][N]), &ONE, &(p[m][2*N]), &ONE);
      cTMP = 1.0 / r[m][2];
      zscal_(&N, &cTMP, &(p[m][2*N]), &ONE);
      cTMP = bnrm*c[m][2]*f[m];
      zaxpy_(&N, &cTMP, &(p[m][2*N]), &ONE, &(x[m][0]), &ONE);
      f[m] = -conj(s[m][2])*f[m];
      h[m] = cabs(-conj(s[m][2]))*h[m];
      if(h[m]/bnrm < 1e-13){
        conv_num++;
        is_conv[m] = j;
      }
      c[m][0]=c[m][1]; c[m][1]=c[m][2];
      s[m][0]=s[m][1]; s[m][1]=s[m][2];
    }
    dTMP = 1.0 / beta[1];
    zdscal_(&N, &dTMP, &(q[2*N]), &ONE);
    beta[0] = beta[1];
    zcopy_(&N, &(q[N]),   &ONE, &(q[0]), &ONE);
    zcopy_(&N, &(q[2*N]), &ONE, &(q[N]), &ONE);

    if(conv_num == M){
      break;
    }
  }
  end_time = time(NULL);
  // Finalize
  fprintf(stdout, "# shifted MINRES method\n");
  fprintf(stdout, "# matrix=%s\n", FNAME);
  fprintf(stdout, "# iteration=%d, status=%d/%d, time=%ld\n", j, conv_num,M, end_time-start_time);
  fprintf(stdout, "# k, real(sigma[k]), imag(sigma[k]) conv_itr REAL_RES TRUE_RES\n");
  for(m=0; m<M; m++){
    SpMV(A_row,A_col,A_ele, &(x[m][0]), res, N);
    cTMP = sigma[m];
    zaxpy_(&N, &cTMP, &(x[m][0]), &ONE, res, &ONE);
    cTMP = -1.0;
    zaxpy_(&N, &cTMP, b, &ONE, res, &ONE);
    fprintf(stdout, "%d %lf %lf %d %e %e\n",
            m, creal(sigma[m]), cimag(sigma[m]), is_conv[m], h[m]/bnrm, dznrm2_(&N,res,&ONE)/bnrm);
  }
  free(A_row);free(A_col);free(A_ele); free(b); free(sigma);
  free(q); free(f); free(h);
  for(m=0; m<M; m++){
    free(r[m]); free(c[m]); free(s[m]);
    free(p[m]); free(x[m]);
  }
  free(r); free(c); free(s); free(p); free(x);
  free(is_conv);
  free(res);
  return 0;
}
