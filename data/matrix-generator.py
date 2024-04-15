import numpy as np
import scipy as sp

N     = 10
fname = "sample-matrix.mtx"
cmnt  = "Generated by matrix-generator.py"

Apre = np.random.randint(-9, 10, (N, N), dtype="float64")
Apst = (Apre + Apre.T) / 2.0

sp.io.mmwrite(fname, A, comment=cmnt, field='real', precision=5, symmetry='symmetric')
