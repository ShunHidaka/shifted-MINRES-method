# shifted-MINRES-method

## Abstract
(standard) Shifted linear systems:
```math
(A + \sigma^{(m)} I) \textbf{x}^{(m)} = \textbf{b},\qquad (m=1,\dots,M),
```
here, $A$ is real symmetric or complex Hermitian.  
This repository is the solver library that offer shifted MINRES method.

# Directory tree
* article
  * Programs that used in "Performance of the shifted minimal residual method for multiply shifted linear systems with real symmetric or complex Hermitian coefficient matrices"
* doc
  * This will the manuals.
  * Underconstruction
* src-c
  * shifted MINRES method in C implementation
  * Please see "src-c/README.md"
* src-matlab
  * shifted MINRES method in MATLAB implementation
  * Please see "src-matlab/README.md"

# Documents
* [Manual for library](https://github.com/ShunHidaka/shifted-MINRES-method/blob/main/doc/manual.pdf) (very old)

# TODO
* update Manual
* update and test matlab
