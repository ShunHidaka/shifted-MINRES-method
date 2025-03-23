# How to use

## Requirement


## How to use
Change *Makefile* ```LAFLAGS``` from ```-lgfortran -lblas -llapack``` to your BLAS/LAPACK.
```bash
make all
```

## Known Issues
- **TRUE_RES always reported as 1.0 (sMINRES)**:
  - Under certain OpenBLAS library versions, the TRUE_RES value in `sminres.out` may incorrectly appear as 1.0.
  - **Cause**: A known bug associated with the `zrotg` function in older versions of OpenBLAS.
    - See: https://github.com/OpenMathLib/OpenBLAS/issues/4909
  - **Solution**:
    - Update OpenBLAS to a newer version that resolves this bug.
    - Use an alternative BLAS library (e.g., Reference-BLAS/LAPACK or Intel MKL).
    - 以下のように sminres.c を書き換えてください

