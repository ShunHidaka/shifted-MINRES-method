SHELL = /bin/bash
CC = gcc
CFLAGS   = -lm -Wall -O3 -fopenmp
LAFLAGS  = -lgfortran -lblas -llapack
ifdef debug
	CFLAGS = -lm -Wall -O0 -g3
endif

.PHONY: all run solver clean

all: run

run: solver
	./sminres.out 1 2> stderr1_sminres.dat
	./scocg.out   1 2> stderr1_scocg.dat
	./sminres.out 2 2> stderr2_sminres.dat
	./scocg.out   2 2> stderr2_scocg.dat
	./sminres.out 3 2> stderr3_sminres.dat
	./sbicg.out   3 2> stderr3_sbicg.dat
	./sminres.out 4 2> stderr4_sminres.dat
	./sbicg.out   4 2> stderr4_sbicg.dat
	gnuplot graph1.gp

solver: sminres.out scocg.out sbicg.out
sminres.out: sminres.c function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)
scocg.out: scocg.c function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)
sbicg.out: sbicg.c function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)

clean:
	rm -f *.out
	rm -f *.dat
