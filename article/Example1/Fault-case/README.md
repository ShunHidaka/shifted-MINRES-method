# Example1: Validation of the setting (selection of optimal seed)

This directory contains programs and input files used in Section 3.1 of the paper  
**"Performance of the shifted minimal residual method for multiply shifted linear systems with real symmetric or complex Hermitian coefficient matrices"**

We investigate the necessity of optimal seed selection in shifted COCG method.  
This experiment demonstrates that inappropriate seeds can lead to solver instability or failure.

---

## Execution
To run all fault-case tests, execute:
```bash
make run
```

Among the solvers used, the scocg(seed: σ=1.0i) implementation outputs residual norms at each iteration.  
The output format is as follows:
```
ITER  ALGO(σ=0.0)  TRUE(σ=0.0)  ALGO(σ=0.01i) TRUE(σ=0.01i) ALGO(σ=0.1i)  TRUE(σ=0.1i)  ALGO(σ=1.0i) TRUE(σ=1.0i)
```
Here:
* "ITER" denotes the number of iterations.
* "ALGO" denotes the residual computed within the algorithm.
* "TRUE" denotes the actual residual $\|\| {\bf b} - (A + \sigma I){\bf x} \|\|$ computed explicitly.

This output helps analyze the convergence behavior of each shifted system and identify failure cases.


## Notes
* All experiments in this directory use the matrix PPE3594, converted to CSR format.  
* See the top-level `README` for details on matrix data preparation.
