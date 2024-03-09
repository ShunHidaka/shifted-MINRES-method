#include "function_sminres.h"

int main(void){
  printf("test\n");
  int N = 10;
  double complex *rhs;
  int M = 10;
  double complex *sigma;
  double threshold = 0.1;
  sminres_init(N, rhs, M, sigma, threshold);
  sminres_finalize();
  return 0;
}
/*
gcc sample2.c function_sminres.c -lm -lgfortran -lblas -llapack
*/
