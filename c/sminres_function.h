#include "definition.h"

void sminres_init(const int input_N, double complex *input_rhs,
                  const int input_M, double complex *input_sigma,
                  double complex *input_v, double complex *x,
                  const double input_threshold);

void sminres_update(int j, double complex *input_v, double complex *input_Av,
                    double complex *x, int *status);

void sminres_finalize();

void sminres_getresidual();
