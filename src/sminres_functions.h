#include "definition.h"

// Initialize shifted MINRES method
void sminres_initialize(const int input_N, double complex *input_rhs,
                        const int input_M, double complex *input_sigma,
                        double complex *input_q,
                        const double input_threshold);

// Update approximate solution
void sminres_update(int j, double complex *input_q, double complex *input_Aq,
                    double complex **x, int *status);

// Finalize shifted MINRES method
void sminres_finalize();

// Return residual norms
void sminres_getresidual(double *input_h);
