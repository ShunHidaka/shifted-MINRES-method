# How to use
## BLAS/LAPACK
Change *Makefile* ```LAFLAGS``` from ```-lgfortran -lblas -llapack``` to your BLAS/LAPACK.

## Prepare matrix file
Matrix data can acquire [ELSES MATRIX LIBRARY](http://www.elses.jp/matrix/).  
Ex.
```bash
wget
tar -xzvf
```
In addition, our implementation use **CSR format**(Compressed Sparse Row).  
This format is superior in terms of memory consumption and allows for high-speed matrix-vector multiplication.  
We are providing the [Python program](https://github.com/ShunHidaka/shifted-MINRES-method/blob/main/article/converter.py) that convert "Matrix Market format" to "CSR format".
Please run:
```python
$ converter.py
file namt: (Matrix Market format file name)
```

## Compile and Run
```bash
make all
```

## Note.
* Segmentation fault in scocg and sbicg
  * check the seed 's' and the number of shifts 'M'
