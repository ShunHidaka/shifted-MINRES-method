# shifted MINRES method

This repository contains the source codes and data used in the paper:  
**Performance of the shifted minimal residual method for multiply shifted linear systems with real symmetric or complex Hermitian coefficient matrices**  
**Authors:** Shuntaro Hidaka, Shuhei Kudo, Takeo Hoshi, Yusaku Yamamoto.

This work implements and evaluates the shifted MINRES, shifted COCG, and shifted BiCG methods for solving multiply shifted linear systems.

---

## Requirements

- C compiler (e.g., `gcc`)
- BLAS and LAPACK libraries (e.g., Reference-BLAS/LAPACK, OpenBLAS, Intel MKL)
- `make`
- python (numpy and scipy)

Tested with:
- gcc (version 9.4.0), Reference-BLAS/LAPACK (version 3.9.0-1build1)
- gcc (version 8.5.0), OpenBLAS (version 0.3.24)
  - OpenBLASを用いる場合は注意が必要です。


## Compilation

Open the `Makefile` and update the following lines according to your environment:
```makefile
CFLAGS  = gcc
LAFLAGS = -lgfortran -lblas -llapack
PYTHON  = python.exe
```
Then download and convert MATRIX data, and compile all programs with:
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
