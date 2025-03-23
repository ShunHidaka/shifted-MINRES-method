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
- python3 with `numpy` and `scipy` (used for matrix data conversion)

Tested with:
- gcc (version 9.4.0) + Reference-BLAS/LAPACK (version 3.9.0-1build1)
- gcc (version 8.5.0) + OpenBLAS (version 0.3.27)
  - **Caution:** when using OpenBLAS prior to version 0.3.27 (see below)

## Directory Structure
This repository consists of the following directories:
- `Example1/`  
  Contains the programs used in Section 3.1 of the paper.  
  See the `/Example1/README.md` for details.
- `Example2/`  
  Contains the programs used in Section 3.2 of the paper.  
  See the `/Example2/README.md` for details.
- `Example3/`  
  Contains the programs used in Section 3.3 of the paper.  
  See the `Example3/README.md`for details.
- `Original/`
  This directory contains the original implementation without any modifications for reproducibility.  
  See the `Original/README.md`for details.
- `MATRIX/`  
  Directory for storing matrices used in numerical experiments. Created automatically by `make`.

## Compilation
Edit the `Makefile` as needed for your environment:
```makefile
CC      = gcc
LAFLAGS = -lgfortran -lblas -llapack
PYTHON  = python.exe
```
Then compile all programs with:
```bash
make solver
```
This compiles all experimental programs in each `Example*/` subdirectory.

## Matrix Data Preparation
Run the following command to download and convert matrix data:
```bash
make init
```
This will download and convert the necessary matrix files into the `MATRIX/` directory.  
If you would like to convert your own matrix file in Matrix Market format to CSR format, use the following command:
```bash
python converter.py input_matrix.mtx output_matrix.csr
```

## Usage
After compilation, move into each example directory and run:
```bash
cd Example1
make run
```
This will execute the experiment under the same conditions as described in the paper.  
Repeat for `Example2/` and `Example3/` as needed.  
For more details about each experiment, please refer to the `README.md` inside each directory.

## Known Issues
- **TRUE_RES always reported as 1.0 (sMINRES)**:
  - Under certain OpenBLAS library versions, the TRUE_RES value in `sminres.out` may incorrectly appear as 1.0.
  - **Cause**: A known bug in the `zrotg` function in OpenBLAS versions prior to 0.3.27.
    - See: https://github.com/OpenMathLib/OpenBLAS/issues/4909
  - **Workarounds**:
    - Update OpenBLAS version 0.3.27 or later.
    - Use an alternative BLAS implementation (e.g., Reference-BLAS/LAPACK or Intel MKL).
    - Optionally, modify the source to use LAPACK's `zlartg` instead of `zrotg` for Givens rotation.  
      Example modification:  
      Replace:
      ```c
      zrotg_(&(r[m][2]), &cTMP, &(c[m][2]), &(s[m][2]));
      ```
      With:
      ```c
      double complex tmp_zlartg = beta[1];
      zlartg_(&(r[m][2]), &cTMP, &(c[m][2]), &(s[m][2]), &tmp_zlartg);
      r[m][2] = tmp_zlartg;
      ```

## Citation
If you use this code, please cite:
``` bibtex
@article{
  author  = {Shuntaro Hidaka, Shuhei Kudo, Takeo Hoshi, Yusaku Yamamoto},
  title   = {Performance of the shifted minimal residual method for multiply shifted linear systems with real symmetric or complex Hermitian coefficient matrices},
  doi     = {}, % to be updated
  journal = {Computer Physics Communications},
  volume  = {}, % to be updated
  pages   = {}, % to be updated
  year    = {} % to be updated
}
```
