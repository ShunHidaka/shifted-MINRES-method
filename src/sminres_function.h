#include "definition.h"

// shifted MINRES法の開始
void sminres_init(const int input_N, double complex *input_rhs,
                  const int input_M, double complex *input_sigma,
                  double complex *input_v, double complex *x,
                  const double input_threshold);

// 近似解の更新
void sminres_update(int j, double complex *input_v, double complex *input_Av,
                    double complex *x, int *status);

// shifted MINRES法の終了
void sminres_finalize();

// Lanczos過程で得られる定数alpha, betaを返す
void sminres_getresidual(double *input_h);
