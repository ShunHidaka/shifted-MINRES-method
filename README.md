# shifted-MINRES-method

## 概要
標準シフト線形方程式
```math
(A + \sigma_k I) \textbf{x} = \textbf{b},\qquad (k=1,\dots,M).
```
に対する shifted MINRES法 を提供するソルバーライブラリ．
ただし， $A$ は実対称・複素エルミート行列とする．


## 使用要件
* Cコンパイラ
* Fortranコンパイラ
* BLAS
* LAPACK

## 使用法
### sminres_solver.c
行列ベクトル積をソルバー側で行う．

### sminres_function.c
行列ベクトル積を使用者側で行う．
