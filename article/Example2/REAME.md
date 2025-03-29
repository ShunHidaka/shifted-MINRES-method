# Example2: Performance evaluation

This directory contains programs and input files used in Section 3.2 of the paper  
**"Performance of the shifted minimal residual method for multiply shifted linear systems with real symmetric or complex Hermitian coefficient matrices"**

We evaluate the performance of shifted Krylov subspace methods for complex Hermitian systems, using complex shifts defined as:

$$
\sigma^{(m)} = 0.01 \exp\left(\frac{2\pi i}{10}(m-0.5)\right), \quad (m=1,\dots,10)
$$

## Execution
To run the experiments, execute:
```bash
make run
```
This will execute performance evaluations and store results in corresponding subdirectories listed below.

## Directory Structure
- **add-Exp/**  
  Experiments with shifts defined as $\sigma^{(m)} = 0.1 \exp\left(\frac{2\pi i}{10}(m-0.5)\right), \quad (m=1,\dots,10)$  
  Used to obtain numerical results reported in Tables 9 and 10 of the paper.
- **residual/**  
  Contains scripts and data used to generate Figure 1, illustrating residual behaviors.
- **time-iter/**  
  Contains performance data used for Tables 6, 7, and 8, evaluating computational time and iteration counts.

## Notes
- All experiments use the complex Hermitian matrix VCNT900h, converted to CSR format.
- For more details on matrix preparation and project structure, see the top-level [README](../README.md).
