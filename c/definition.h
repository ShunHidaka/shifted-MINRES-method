#ifndef INCLUDE_DEFINITION_H
#define INCLUDE_DEFINITION_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <cblas.h>

// x <= a*x
void zdscal_(int *n, double *a, double complex *x, int *incx);
void zscal_(int *n, double complex *a, double complex *x, int *incx);

// y <= alpha*x + y
void zaxpy_(int *n, double complex *alpha, double complex *x, int *incx, double complex *y, int *incy);

// y <= x
void zcopy_(int *n, double complex *x, int *incx, double complex *y, int *incy);

// ret <= <x, y> = x^H y
double complex zdotc_(int *n, double complex *x, int *incx, double complex *y, int *incy);
// ret <= <x, y>' = x^T y
double complex zdotu_(int *n, double complex *x, int *incx, double complex *y, int *incy);

// ret <= ||x||_2
double dznrm2_(int *n, double complex *x, int *incx);

// 倍精度複素数で与えられた点に対してy座標をゼロにするGivens回転を求める
void zrotg_(double complex *a, double complex *b, double *c, double complex *s);

// 倍精度複素数で生成されたGivens回転を適用する
void zrot_(int *n, double complex *x, int *incx, double complex *y, int *incy);

/*
typedef struct{
  double r;
  double i;
} doublecomplex;

// x <= a*x
void zdscal_(int *n, double *a, doublecomplex *x, int *incx);
void zscal_(int *n, doublecomplex *a, doublecomplex *x, int *incx);

// y <= alpha*x + y
void zaxpy_(int *n, doublecomplex *alpha, doublecomplex *x, int *incx, doublecomplex *y, int *incy);

// y <= x
void zcopy_(int *n, doublecomplex *x, int *incx, doublecomplex *y, int *incy);

// ret <= <x, y> = x^H y
doublecomplex zdotc_(int *n, doublecomplex *x, int *incx, doublecomplex *y, int *incy);
// ret <= <x, y>' = x^T y
doublecomplex zdotu_(int *n, doublecomplex *x, int *incx, doublecomplex *y, int *incy);

// ret <= ||x||_2
double dznrm2_(int *n, doublecomplex *x, int *incx);

// 倍精度複素数で与えられた点に対してy座標をゼロにするGivens回転を求める
void zrotg_(doublecomplex *a, doublecomplex *b, double *c, doublecomplex *s);

// 倍精度複素数で生成されたGivens回転を適用する
void zrot_(int *n, doublecomplex *x, int *incx, doublecomplex *y, int *incy);
*/
#endif
