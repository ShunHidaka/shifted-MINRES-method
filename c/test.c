#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <cblas.h>
#include "definition.h"

int main(void){
  printf("TEST\n");
  doublecomplex a = {1.0, 2.0};
  doublecomplex b = a;
  doublecomplex c;
  c.r = a.r; c.i = a.i;
  printf("%lf %lf\n", a.r, a.i);
  printf("%lf %lf\n", b.r, b.i);
  printf("%lf %lf\n", c.r, c.i);
  printf("%p %p %p\n", &a, &b, &c);

  double complex d = 1.0+2.0I;
  double complex e = 1.0+2.0I;
  printf("%p %p %p\n", &d, &e, &d);

  void *temp1 = &d;
  double *temp2 = temp1;

  printf("%lf %lf\n", temp2[0], temp2[1]);
  printf("%p %p\n", &(temp2[0]), &(temp2[1]));
  return 0;
}
