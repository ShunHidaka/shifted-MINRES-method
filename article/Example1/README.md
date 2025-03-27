# Example1: Performance Evaluation with Different Right-Hand Side Vectors

This directory contains programs and input files used in Section 3.1 of the paper  
**"Performance of the shifted minimal residual method for multiply shifted linear systems..."**

We evaluate the performance of the shifted MINRES method by varying the right-hand side (RHS) vector while solving multiply shifted linear systems with the **PPE3594** matrix (real symmetric).

---

## Execution
To run all tests in this example, execute:
```bash
make run
```
This will perform the shifted MINRES computation for three types of RHS vectors:

| RHS Type          | Description                       | Input Argument |
|-------------------|-----------------------------------|----------------|
| All-one vector	| A vector where all entries are 1	| 0              |
| Random vector 1	| `srand(1)` seeded randomness	    | 1              |
| Random vector 2	| `srand(2)` seeded randomness	    | 2              |


## Notes
* All experiments in this directory use the matrix PPE3594, converted to CSR format.  
* See the top-level `README` for details on matrix data preparation.
