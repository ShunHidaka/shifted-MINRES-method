#include "definition.h"
#include "function_util.h"
char *FNAME = "/home/stern/multi/data/VCNT900h.csr";

int main(int argc, char *argv)
{
  int i;
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

  /* set options and allocate memory */
  itermax   = N;
  threshold = 1e-12;
  status    = (int *)calloc(3 , sizeof(int));
  dsfmt_init_gen_rand(&dsfmt, 120);
  
  /* generate points for integration and rhs vector */
  nz    = 10; // 積分点の数
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

  /* shifted MINRES method */
  x = (double complex *)calloc(nz*N, sizeof(double complex));
  solve_sminres(z, A_row,A_col,A_ele, rhs, x, nz, N, itermax,threshold,status);

  {
    double res;
    double complex *tmp = (double complex *)calloc(N, sizeof(double complex));
    for(int k=0; k<nz; k++){
      SpMV(A_row,A_col,A_ele, &(x[k*N]), tmp, N);
      for(i=0; i<N; i++) tmp[i] = rhs[i] - (tmp[i] + z[k]*x[k*N+i]);
      int ONE = 1;
      res = dznrm2_(&N, tmp, &ONE);
      fprintf(stdout, "%d %e\n", k, res);
    }
    free(tmp);
  }

  return 0;
}
