#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#undef I
#define Im_I _Complex_I
#include "mmio.h"

// Read Matrix Market format & Convert BLAS Packed Storage
// only "matrix cordinate complex hermitian"
// or   "matrix cordinate real symmetric"
void read_mtx(const char *fname, int *N, double complex **matrix);

// Set shifts
void set_shifts(int *M, double complex **shifts);

// Set Right-Hand-Side vector
void set_rhs(const int N, double complex **rhs);
