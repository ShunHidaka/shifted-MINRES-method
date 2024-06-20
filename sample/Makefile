CC = gcc
CFLAGS  = -lm -Wall -Wextra -O2
#LAFLAGS = -lgfortran -lblas -llapack
LAFLAGS = -I/opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib -lopenblas

all: sample1.out sample2.out sample3.out
sample1.out: sample1.c functions.c ../src/c/sminres_solver.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)
sample2.out: sample2.c functions.c ../src/c/sminres_function.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)

sample3.out: sample3.c functions.c ../src/c/scocg_solver.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)

clean:
	rm *.out
