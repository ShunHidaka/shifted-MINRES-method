import numpy as np
import scipy as sp

matrix = sp.io.mmread("ELSES_MATRIX_PPE3594_20160426_A.mtx").tocsr()
shift = sp.sparse.csr_matrix(0.1 * 1j * np.identity(matrix.shape[0], dtype=np.complex128))
A = matrix + shift

N = len(matrix.indptr)-1
b = np.empty(N, np.complex128)
for i in range(N):
    b[i] = 1.0

x = np.empty(N, np.complex128)
for i in range(10):
    x = sp.sparse.linalg.spsolve(A, b)

fp = open("answer.dat", 'w', encoding="utf-8")
for i in range(N):
    print("%lf %lf\n" % (x[i].real, x[i].imag), end=" ", file=fp)
