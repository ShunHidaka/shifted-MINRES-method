# shifted-MINRES-method

## Abstract
(standard) Shifted linear systems:
```math
(A + \sigma^{(m)} I) \textbf{x}^{(m)} = \textbf{b},\qquad (m=1,\dots,M),
```
here, $A$ is real symmetric or complex Hermitian.  
This repository is the solver library that offer shifted MINRES method.


## Requirement
* C compiler (tested only GCC)
* BLAS (tested Reference-BLAS and Open BLAS)
* LAPACK

Note: This software is developed and tested WSL2(Ubuntu)

## How to use
This softere provided two-type solver
### sminres_function.c
* sminres_initalize()
* sminres_update()
* sminres_finalize()
* sminres_getresidual()
