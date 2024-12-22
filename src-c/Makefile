CC = gcc
CFLAGS  = -lm -Wall -Wextra -O2
LAFLAGS = -lgfortran -lblas -llapack
#LAFLAGS = -I/opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib -lopenblas

all: sample.out
sample.out: sample.c function_util.c sminres_functions.c mmio.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)

clean:
	rm *.out
