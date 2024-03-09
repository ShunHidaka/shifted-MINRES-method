#include "definition.h"
#include "function_sminres.h"

int N, M;
double threshold;
double r0_nrm;
double complex *v;
double *alpha, *beta;
double complex *T;
double *G1; double complex *G2;
double complex *p, *f; double *h;
void sminres_init(const int input_N, double complex *input_rhs,
                  const int input_M, double complex *input_sigma,
                  const double input_threshold)
{
  N         = input_N;
  M         = input_M;
  threshold = input_threshold;
  v      = (double complex *)calloc(N*3,   sizeof(double complex));
  alpha  = (double         *)calloc(1,     sizeof(double));
  beta   = (double         *)calloc(2,     sizeof(double));
  T      = (double complex *)calloc(M*4,   sizeof(double complex));
  G1     = (double         *)calloc(M*3,   sizeof(double));
  G2     = (double complex *)calloc(M*3,   sizeof(double complex));
  p      = (double complex *)calloc(M*N*3, sizeof(double complex));
  f      = (double complex *)calloc(M,     sizeof(double complex));
  h      = (double         *)calloc(M,     sizeof(double));
  
}
void sminres_update()
{
  
}
void sminres_finalize()
{
  free(v);
  free(alpha); free(beta);
  free(T); free(G1); free(G2);
  free(p); free(f); free(h);
}
