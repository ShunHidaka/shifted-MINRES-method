# shifted MINRES method in MATLAB
Solver library for the standard shifted linear systems:
```math
(A + \sigma^{(m)} I) \textbf{x}^{(m)} = \textbf{b},\qquad (m=1,\dots,M),
```

*matrix A is Real Symmetric*

## Requirement
* MATLAB
* GNU Octave (Not tested, but maybe runnable)

## How to use
```Matlab
function [x, flag, rres, itrs] = shifted_minres(A, rhs, N, sigma, M, max_itr, threshold)
% Argument：
%   "Matrix" A, "Right-Hand-Side-vector" rhs, "Matrix-size" N,
% 　"Shifts" sigma, "Number-of-shits" M,
%   "Maximum-iteration" max_itr, "relative-residual-norm-threshold" threshold
% Return：
%   "Approximate-solution" x (M \times N)
%   "Status" flag
%   "Relative-residual-norm" rres (M)
%   "Converged iterations" itrs (M)
```

## TODO
* Implement the case of $A$ is complex Hermitian
* Enough debug and test !!
