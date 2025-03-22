SHELL = /bin/bash
CC = gcc
CFLAGS   = -lm -Wall -O3 -fopenmp
LAFLAGS  = -lgfortran -lblas -llapack
#LAFLAGS  = -I/opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib -lopenblas
PYTHON   = python.exe
ifdef debug
	CFLAGS = -lm -Wall -O0 -g3
endif

.PHONY: all init solver example1 Ex1 Fau add res tim Ex3 clean allclean

all: init solver

Ex1 = Example1
Fau = $(Ex1)/Fault-case
Ex2 = Example2
add = $(Ex2)/add-Exp
res = $(Ex2)/residual
tim = $(Ex2)/time-iter
Ex3 = Example3
sm = sminres.out
sc = scocg.out
sb = sbicg.out

##########
solver: sample example1 example2 example3

sample: $(sm) $(sc) $(sb)
$(sm): sminres.c function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)
$(sc): scocg.c function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)
$(sb): sbicg.c function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)

example1: Ex1 Fau
Ex1: $(Ex1)/$(sm) $(Ex1)/$(sc)
$(Ex1)/$(sm): $(Ex1)/sminres.c $(Ex1)/function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)
$(Ex1)/$(sc): $(Ex1)/scocg.c   $(Ex1)/function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)
Fau: $(Fau)/$(sm) $(Fau)/$(sc)
$(Fau)/$(sm): $(Fau)/sminres.c $(Fau)/function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)
$(Fau)/$(sc): $(Fau)/scocg.c   $(Fau)/function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)

example2: add res tim
add: $(add)/$(sm) $(add)/$(sc) $(add)/$(sb)
$(add)/$(sm): $(add)/sminres.c $(add)/function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)
$(add)/$(sc): $(add)/scocg.c $(add)/function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)
$(add)/$(sb): $(add)/sbicg.c $(add)/function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)
res: $(res)/$(sm) $(res)/$(sc) $(res)/$(sb)
$(res)/$(sm): $(res)/sminres.c $(res)/function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)
$(res)/$(sc): $(res)/scocg.c $(res)/function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)
$(res)/$(sb): $(res)/sbicg.c $(res)/function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)
tim: $(tim)/$(sm) $(tim)/$(sc) $(tim)/$(sb)
$(tim)/$(sm): $(tim)/sminres.c $(tim)/function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)
$(tim)/$(sc): $(tim)/scocg.c $(tim)/function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)
$(tim)/$(sb): $(tim)/sbicg.c $(tim)/function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)

example3: $(Ex3)/$(sm) $(Ex3)/$(sc) $(Ex3)/$(sb)
$(Ex3)/$(sm): $(Ex3)/sminres.c $(Ex3)/function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)
$(Ex3)/$(sc): $(Ex3)/scocg.c $(Ex3)/function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)
$(Ex3)/$(sb): $(Ex3)/sbicg.c $(Ex3)/function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)

##########
init: MATRIX/PPE3594_A.csr MATRIX/CLIQ6912std_A.csr MATRIX/CLIQ55296std_A.csr MATRIX/VCNT900h_A.csr MATRIX/VCNT10800h_A.csr
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
	wget http://www.damp.tottori-u.ac.jp/~hoshi/elses_matrix/ELSES_MATRIX_VCNT900h_20130501.tgz
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

##########
clean:
	rm -f *.out $(Ex1)/*.out $(Fau)/*.out $(add)/*.out $(res)/*.out $(tim)/*.out $(Ex3)/*.out
	rm -f *.dat $(Ex1)/*.dat $(Fau)/*.dat $(add)/*.dat $(res)/*.dat $(tim)/*.dat $(Ex3)/*.dat
allclean: clean
	rm -f MATRIX/*
	rmdir MATRIX/
