SHELL = /bin/bash
CC = gcc
CFLAGS   = -lm -Wall -O3 -fopenmp
LAFLAGS  = -lgfortran -lblas -llapack
ifdef debug
	CFLAGS = -lm -Wall -O0 -g3
endif

.PHONY: all init clean

all: shift

shift: sminres.out scocg.out sbicg.out
sminres.out: sminres.c function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)
scocg.out: scocg.c function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)
sbicg.out: sbicg.c function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)

init:
	mkdir MATRIX
	wget -P MATRIX http://www.damp.tottori-u.ac.jp/~hoshi/elses_matrix/ELSES_MATRIX_PPE3594_20160426.tgz
	tar -xzvf ./MATRIX/ELSES_MATRIX_PPE3594_20160426.tgz
	python converter.py ./MATRIX/ELSES_MATRIX_PPE3594_20160426/ELSES_MATRIX_PPE3594_A ./MATRIX/PPE3594_A

clean:
	rm -f *.out *.dat
allclean: clean
	rm -rf MATRIX
