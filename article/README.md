# shifted MINRES method

This repository contains the source codes and data used in the paper:  
**Performance of the shifted minimal residual method for multiply shifted linear systems with real symmetric or complex Hermitian coefficient matrices**  
**Authors:** Shuntaro Hidaka, Shuhei Kudo, Takeo Hoshi, Yusaku Yamamoto.

The objective of this study is to propose the shifted MINRES method as a solver for multiply shifted linear systems with real symmetric or complex Hermitian coefficient matrices, and to evaluate its performance.

---

## Requirements
- C compiler (e.g., `gcc`)
- BLAS and LAPACK libraries (e.g., Reference-BLAS/LAPACK, OpenBLAS, Intel MKL)
- `make`
- python3 with `numpy` and `scipy` (used for matrix data conversion)

Tested with:
- gcc (version 9.4.0) + Reference-BLAS/LAPACK (version 3.9.0-1build1)
- gcc (version 8.5.0) + OpenBLAS (version 0.3.27)
  - **Caution:** when using OpenBLAS prior to version 0.3.27 ([see below](https://github.com/ShunHidaka/shifted-MINRES-method/tree/main/article#known-issues))

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

### Test Matrices
The following matrices from the [ELSES matrix library](http://www.elses.jp/matrix/) are used in our experiments:

| Matrix Name      | Size    | NNZ        | Type              | Link                                     |
|------------------|---------|------------|-------------------|------------------------------------------|
| PPE3594          | 3594    | 79,968     | Real symmetric    | http://www.elses.jp/matrix/#PPE3594      |
| CLIQ6912std      | 6912    | 208,544    | Real symmetric    | http://www.elses.jp/matrix/#CLIQ6912std  |
| CLIQ55296std     | 55296   | 1,652,352  | Real symmetric    | http://www.elses.jp/matrix/#CLIQ55296std |
| VCNT900h         | 900     | 683,892    | Complex Hermitian | http://www.elses.jp/matrix/#VCNT900h     |
| VCNT10800h       | 10800   | 8,511,588  | Complex Hermitian | http://www.elses.jp/matrix/#VCNT10800h   |

These matrices are provided in [Matrix Market format](https://math.nist.gov/MatrixMarket/), a widely used format for sparse matrices.  
To enable **efficient matrix-vector multiplication**, the matrices are **converted to CSR (Compressed Sparse Row) format** using the Python script `converter.py`.
The converted text-based `.csr` files are stored in the `MATRIX/` directory and are read directly by the solver programs.

If, for any reason, the `make init` command does not work properly, you can manually download and convert the matrix files by running the following commands:
```bash
# Make the MATRIX/ directory
mkdir MATRIX
# Download matrix archives from the ELSES matrix library
wget http://www.damp.tottori-u.ac.jp/~hoshi/elses_matrix/ELSES_MATRIX_PPE3594_20160426.tgz
wget http://www.damp.tottori-u.ac.jp/~hoshi/elses_matrix/ELSES_MATRIX_CLIQ6912std_20130109.tgz
wget http://www.damp.tottori-u.ac.jp/~hoshi/elses_matrix/ELSES_MATRIX_CLIQ55296std_20130109.tgz
wget http://www.damp.tottori-u.ac.jp/~hoshi/elses_matrix/ELSES_MATRIX_VCNT900h_20130501.tgz
wget http://www.damp.tottori-u.ac.jp/~hoshi/elses_matrix/ELSES_MATRIX_VCNT10800h_20130501.tgz
# Extract each archive
tar -xzvf ./ELSES_MATRIX_PPE3594_20160426.tgz
tar -xzvf ./ELSES_MATRIX_CLIQ6912std_20130109.tgz
tar -xzvf ./ELSES_MATRIX_CLIQ55296std_20130109.tgz
tar -xzvf ./ELSES_MATRIX_VCNT900h_20130501.tgz
tar -xzvf ./ELSES_MATRIX_VCNT10800h_20130501.tgz
# Convert Matrix Market format to CSR format using the provided Python script
python converter.py ./ELSES_MATRIX_PPE3594_20160426/ELSES_MATRIX_PPE3594_20160426_A.mtx MATRIX/PPE3594_A.csr
python converter.py ./ELSES_MATRIX_CLIQ6912std_20130109/ELSES_MATRIX_CLIQ6912std_A.mtx MATRIX/CLIQ6912std_A.csr
python converter.py ./ELSES_MATRIX_CLIQ55296std_20130109/ELSES_MATRIX_CLIQ55296std_A.mtx MATRIX/CLIQ55296std_A.csr
python converter.py ./ELSES_MATRIX_VCNT900h_20130501/ELSES_MATRIX_VCNT900h_A.mtx MATRIX/VCNT900h_A.csr
python converter.py ./ELSES_MATRIX_VCNT10800h_20130501/ELSES_MATRIX_VCNT10800h_A.mtx MATRIX/VCNT10800h_A.csr
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
@article{TEMPORARY
  author  = {Shuntaro Hidaka, Shuhei Kudo, Takeo Hoshi, Yusaku Yamamoto},
  title   = {Performance of the shifted minimal residual method for multiply shifted linear systems with real symmetric or complex Hermitian coefficient matrices},
  doi     = {https://doi.org/10.1016/j.cpc.2025.109679},
  journal = {Computer Physics Communications},
  volume  = {}, % to be updated
  pages   = {109679},
  year    = {2025}
}
```
