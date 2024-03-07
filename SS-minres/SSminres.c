#include "definition.h"
#include "function_util.h"
char *FNAME = "/home/stern/multi/data/VCNT900h.csr";

int main(int argc, char *argv)
{
  int i, j, k;
  int nz;
  double gamma, rho;
  double complex *z, *x, *rhs;
  int itermax, *status;
  double threshold;
  dsfmt_t dsfmt;

  /* read R-Symmetric of C-Hermitian Matrix 'A' */
  int N, DATASIZE;
  int *A_row, *A_col;
  double complex *A_ele;
  read_csr(FNAME, &N, &DATASIZE, &A_row, &A_col, &A_ele);
  for(int i=0; i<DATASIZE; i++){
    A_ele[i] = -A_ele[i];
  }
  
  /* generate points for integration and rhs vector */
  dsfmt_init_gen_rand(&dsfmt, 120);
  nz    = 100; // 積分点の数
  gamma = 1.0;
  rho   = 0.2;
  z     = (double complex *)calloc(nz, sizeof(double complex));
  for(i=0; i<nz; i++){
    z[i] = gamma + rho*cexp( 2*PI*I*(i+0.5)/nz );
  }
  rhs = (double complex *)calloc(N, sizeof(double complex));
  for(i=0; i<N; i++){
    rhs[i]  = 2*(dsfmt_genrand_close_open(&dsfmt) - 0.5);
    rhs[i] += 2*(dsfmt_genrand_close_open(&dsfmt) - 0.5)*I;
  }

  /* set options and shifted MINRSE method */
  itermax   = N;
  threshold = 1e-12;
  status    = (int *)calloc(3 , sizeof(int));
  x         = (double complex *)calloc(nz*N, sizeof(double complex));
  solve_sminres(z, A_row,A_col,A_ele, rhs, x, nz, N, itermax,threshold,status);

  /* Numercal Integration in complex plane */
  double complex **s, *tmp_s;
  for(j=0; j<____; j++){
    for(i=0; i<N; i++){
      tmp = 0.0;
      for(k=0; k<nz; k++)
        tmp += x[k*N+i] * cpow( (z[k]-gamma)/rho, 1.0*j+1.0 );
      tmp = tmp / nz;
      S[i][j + ]
    }
  }

  /* Singular Value Decomposition and determine rank 's' */
  ZSVD(N, nk*nr);

  /* projection 'A' to tilde{A} */

  /* diagnolize tilde{A} */

  /* output Eigenvalues of tilde{A} = 'A' */

  free(x); free(rhs); free(z);
  free(status);
  free(A_row); free(A_col); free(A_ele);
  return 0;
}
