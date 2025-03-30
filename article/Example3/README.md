# Example3: Performance evaluation

This directory contains programs and input files used in Section 3.3 of the paper  
**"Performance of the shifted minimal residual method for multiply shifted linear systems with real symmetric or complex Hermitian coefficient matrices"**

We evaluate the performance of shifted Krylov subspace methods for complex Hermitian systems using complex shifts defined as:

$$
\sigma^{(m)} = (-0.501 + 0.001m) + 0.001i, \quad (m = 1, \dots, M)
$$

## Execution
To run the experiments, execute:
```bash
make run
```



## Notes
- All experiments use the complex Hermitian matrix VCNT10800h, converted to CSR format.
- For more details on matrix preparation and project structure, see the top-level [README](../../README.md).
