#ifndef INCLUDE_DEFINITION_H
#define INCLUDE_DEFINITION_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

// x <= a*x
void zdscal_(int *n, double *a, double complex *x, int *incx);
void zscal_(int *n, double complex *a, double complex *x, int *incx);

// y <= alpha*x + y
void zaxpy_(int *n, double complex *alpha, double complex *x, int *incx, double complex *y, int *incy);

// y <= x
void dcopy_(int *n, double *x, int *incx, double *y, int *incy);
void zcopy_(int *n, double complex *x, int *incx, double complex *y, int *incy);

// ret <= <x, y> = x^H y
double complex zdotc_(int *n, double complex *x, int *incx, double complex *y, int *incy);
// ret <= <x, y>' = x^T y
double complex zdotu_(int *n, double complex *x, int *incx, double complex *y, int *incy);

// ret <= ||x||_2
double dznrm2_(int *n, double complex *x, int *incx);

// calculate Givens rotation matrix
void zrotg_(double complex *a, double complex *b, double *c, double complex *s);

// apply Givens rotation (double complex)
void zrot_(int *n, double complex *x, int *incx, double complex *y, int *incy, double *c, double complex *s);

#endif
