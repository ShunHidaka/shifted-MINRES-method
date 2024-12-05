# shifted-MINRES-method

## Abstract
(standard) Shifted linear systems:
```math
(A + \sigma_m I) \textbf{x}^{(m)} = \textbf{b},\qquad (m=1,\dots,M),
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
### sminres_solver.c
行列ベクトル積をソルバー側で行う．

### sminres_function.c
行列ベクトル積を使用者側で行う．
