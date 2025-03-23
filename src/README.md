# shifted-MINRES-method

## Abstract
(standard) Shifted linear systems:
```math
(A + \sigma^{(m)} I) \textbf{x}^{(m)} = \textbf{b},\qquad (m=1,\dots,M),
```
here, $A$ is **real symmetric** or **complex Hermitian**.  
This repository is the solver library that offer shifted MINRES method.


## Requirement
* C compiler (tested only GCC)
* BLAS
* LAPACK (tested Reference-LAPACK and Open BLAS)
Note: This software is developed and tested WSL2(Ubuntu)

If you compile and run *sample.c*, require [Matrix Market IO C routines](https://math.nist.gov/MatrixMarket/mmio-c.html).
* mmio.c
* mmio.h
```bash:downloads file
wget http://math.nist.gov/MatrixMarket/mmio/c/mmio.h
wget http://math.nist.gov/MatrixMarket/mmio/c/mmio.c
```


## How to use
```C:sminres_initialize
Initialize 
```
```C:sminres_update
a
```
```C:sminres_finalize
c
```
```C:sminres_getresidual
D
```
