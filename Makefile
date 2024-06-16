

.o: 
	gcc src/c/sminres_function.c src/c/sminres_solver.c -c -lm -O2 -lgfortran -lblas -llapack
	ar rcs libsminres.a sminres_function.o sminres_solver.o
