CC = gcc
CFLAGS  = -lm -Wall -Wextra -O3
LAFLAGS = -lgfortran -lblas -llapack
#LAFLAGS = -I/opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib -lopenblas

all: sample.out
sample.out: sample.c function_util.c sminres_functions.c mmio.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAGS)

libsminres.a: sminres_functions.c
	$(CC) $^ -c -lm $(CFLAGS) $(LAFLAGS)
	ar rcs libsminres.a sminres_functions.o
	rm sminres_functions.o
sample: sample.c function_util.c mmio.c
	$(CC) $^ -o $@  -L. -lsminres $(CFLAGS) $(LAFLAGS)

clean:
	rm sample.out libsminres.a
