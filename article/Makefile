SHELL = /bin/bash
CC = gcc
CFLAGS   = -lm -Wall -O3 -fopenmp
LAFLAGS  = -lgfortran -lblas -llapack
PYTHON   = python.exe
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

init: MATRIX/PPE3594_A.csr MATRIX/CLIQ6912std_A.csr MATRIX/CLIQ55296std_A.csr
MATRIX/PPE3594_A.csr: | MATRIX
	wget http://www.damp.tottori-u.ac.jp/~hoshi/elses_matrix/ELSES_MATRIX_PPE3594_20160426.tgz
	tar -xzvf ./ELSES_MATRIX_PPE3594_20160426.tgz
	$(PYTHON) converter.py ./ELSES_MATRIX_PPE3594_20160426/ELSES_MATRIX_PPE3594_20160426_A.mtx $@
	rm -rf ELSES_MATRIX_PPE3594_20160426.tgz ELSES_MATRIX_PPE3594_20160426
MATRIX/CLIQ6912std_A.csr: | MATRIX
	wget http://www.damp.tottori-u.ac.jp/~hoshi/elses_matrix/ELSES_MATRIX_CLIQ6912std_20130109.tgz
	tar -xzvf ./ELSES_MATRIX_CLIQ6912std_20130109.tgz
	$(PYTHON) converter.py ./ELSES_MATRIX_CLIQ6912std_20130109/ELSES_MATRIX_CLIQ6912std_A.mtx $@
	rm -rf ELSES_MATRIX_CLIQ6912std_20130109.tgz ELSES_MATRIX_CLIQ6912std_20130109
MATRIX/CLIQ55296std_A.csr: | MATRIX
	wget http://www.damp.tottori-u.ac.jp/~hoshi/elses_matrix/ELSES_MATRIX_CLIQ55296std_20130109.tgz
	tar -xzvf ./ELSES_MATRIX_CLIQ55296std_20130109.tgz
	$(PYTHON) converter.py ./ELSES_MATRIX_CLIQ55296std_20130109/ELSES_MATRIX_CLIQ55296std_A.mtx $@
	rm -rf ELSES_MATRIX_CLIQ55296std_20130109.tgz ELSES_MATRIX_CLIQ55296std_20130109
MATRIX/VCNT900h_A.csr: | MATRIX
	http://www.damp.tottori-u.ac.jp/~hoshi/elses_matrix/ELSES_MATRIX_VCNT900h_20130501.tgz
	tar -xzvf ./ELSES_MATRIX_VCNT900h_20130501.tgz
	$(PYTHON) converter.py ./ELSES_MATRIX_VCNT900h_20130501/ELSES_MATRIX_VCNT900h_A.mtx $@
	rm -rf ELSES_MATRIX_VCNT900h_20130501.tgz ELSES_MATRIX_VCNT900h_20130501
MATRIX/VCNT10800h_A.csr: | MATRIX
	wget http://www.damp.tottori-u.ac.jp/~hoshi/elses_matrix/ELSES_MATRIX_VCNT10800h_20130501.tgz
	tar -xzvf ./ELSES_MATRIX_VCNT10800h_20130501.tgz
	$(PYTHON) converter.py ./ELSES_MATRIX_VCNT10800h_20130501/ELSES_MATRIX_VCNT10800h_A.mtx $@
	rm -rf ELSES_MATRIX_VCNT10800h_20130501.tgz ELSES_MATRIX_VCNT10800h_20130501
MATRIX:
	mkdir -p MATRIX

clean:
	rm -f *.out *.dat
allclean: clean
	rm -rf MATRIX
