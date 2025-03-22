SHELL = /bin/bash
CC = gcc
CFLAGS   = -lm -Wall -O3 -fopenmp
LAFLAGS  = -lgfortran -lblas -llapack
ifdef debug
	CFLAGS = -lm -Wall -O0 -g3
endif

.PHONY: all run clean

all: run

run: shift
	./sminres.out 1 > stdout1_sminres.dat
	./scocg.out   1 > stdout1_scocg.dat
	./sminres.out 2 > stdout2_sminres.dat
	./scocg.out   2 > stdout2_scocg.dat
	./sminres.out 3 > stdout3_sminres.dat
	./sbicg.out   3 > stdout3_sbicg.dat
	./sminres.out 4 > stdout4_sminres.dat
	./sbicg.out   4 > stdout4_sbicg.dat
	gnuplot graph2.gp

shift: sminres.out scocg.out sbicg.out
sminres.out: sminres.c function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)
scocg.out: scocg.c function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)
sbicg.out: sbicg.c function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)

clean:
	rm -f *.out
	rm -f *.dat
