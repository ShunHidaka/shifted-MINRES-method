CC = gcc
CFLAGS  = -lm -Wall -Wextra -O2
LAFLAGS = -lgfortran -lblas -llapack

all: sample1.out sample2.out
sample1.out: sample1.c functions.c ../src/c/sminres_solver.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)

sample2.out: sample2.c functions.c ../src/c/sminres_function.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)
clean:
	rm sample1.out sample2.out
