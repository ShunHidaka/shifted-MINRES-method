SHELL = /bin/bash
CC    = gcc
DSFMT = -DDSFMT_MEXP=19937 -DHAVE_SSE2 dSFMT.c
MATH  = -lm -lgfortran -lblas -llapack -g

a.out: SSminres.c function_util.c
	$(CC) $^ $(DSFMT) -o $@ -O2 $(MATH)
test.out: test.c function_util.c
	$(CC) $^ $(DSFMT) -o $@ $(MATH)
