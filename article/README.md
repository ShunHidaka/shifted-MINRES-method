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
- python3 with `numpy` and `scipy`
  - 行列データを変換するのに使用

Tested with:
- gcc (version 9.4.0) + Reference-BLAS/LAPACK (version 3.9.0-1build1)
- gcc (version 8.5.0) + OpenBLAS (version 0.3.24)
  - **Caution:** when using OpenBLAS prior to version 0.3.27 (see below)

## Directory Structure
This repository consists of the following directories:
- `Example1/`  
  3.1節で使用したプログラムがある
  詳細は そのディレクトリ内の README を
- `Example2/`  
  3.2節で使用したプログラムがある
  詳細は そのディレクトリ内の README を
- `Example3/`  
  3.3節で使用したプログラムがある
  詳細は そのディレクトリ内の README を
- `MATRIX/`  
  数値実験で使用する行列を保存するディレクトリ
  makeで生成される

## Compilation
Open the `Makefile` and update the following variables according to your environment:
```makefile
CC      = gcc
LAFLAGS = -lgfortran -lblas -llapack
PYTHON  = python.exe
```
Then compile all programs with:
```bash
make solver
```

## Matrix Data Preparation
```bash
make init
```

## Known Issues
- **TRUE_RES always reported as 1.0 (sMINRES)**:
  - Under certain OpenBLAS library versions, the TRUE_RES value in `sminres.out` may incorrectly appear as 1.0.
  - **Cause**: A known bug in the `zrotg` function in OpenBLAS versions prior to 0.3.27.
    - See: https://github.com/OpenMathLib/OpenBLAS/issues/4909
  - **Workarounds**:
    - Update OpenBLAS version 0.3.27 or later.
    - Use an alternative BLAS implementation (e.g., Reference-BLAS/LAPACK or Intel MKL).
    - Optionally, modify the source to use LAPACK's `zlartg` instead of `zrotg` for Givens rotation.
