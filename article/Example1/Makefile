SHELL = /bin/bash
CC = gcc
CFLAGS   = -lm -Wall -O3 -fopenmp
LAFLAGS  = -lgfortran -lblas -llapack
ifdef debug
	CFLAGS = -lm -Wall -O0 -g3
endif

.PHONY: all run shift clean

all: run
run: shift
	./sminres.out 0 > rhs0_sminres.dat
	./scocg.out   0 > rhs0_scocg.dat
	./sbicg.out   0 > rhs0_sbicg.dat
	./sminres.out 1 > rhs1_sminres.dat
	./scocg.out   1 > rhs1_scocg.dat
	./sbicg.out   1 > rhs1_sbicg.dat
	./sminres.out 2 > rhs2_sminres.dat
	./scocg.out   2 > rhs2_scocg.dat
	./sbicg.out   2 > rhs2_sbicg.dat

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
