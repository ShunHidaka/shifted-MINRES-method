# Krylov subspace methods for shifted linear systems
Solver library for the standard shifted linear systems:
```math
(A + \sigma^{(m)} I) \textbf{x}^{(m)} = \textbf{b},\qquad (m=1,\dots,M),
```
* shifted CG method  
  for real symmetric / complex Hermitian **positive definite** matrix  + real shifts
* shifted MINRES method  
  for real symmetric / complex Hermitian matrix + real / complex shifts


## TODO
* shifted CG法にseed switchingを実装する
* ドキュメントを作成する
