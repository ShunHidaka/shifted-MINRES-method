#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
/* Resolving variable definition conflicts */
#undef I
#define Im_I _Complex_I
#include "mmio.h"

int main(void){
  int ret_code;
  MM_typecode matcode;
  FILE *f;
  int M, N, NNZ;
  double complex *matrix;

  if((f=fopen("ELSES_MATRIX_DIAB18h_A.mtx", "r")) == NULL){
    fprintf(stderr, "Error: Could not found file\n");
    exit(1);
  }
  if(mm_read_banner(f, &matcode) != 0){
    fprintf(stderr, "Error: Could not process MM banner.\n");
    exit(1);
  }
  if( (mm_is_matrix(matcode)) && (mm_is_sparse(matcode)) &&
      ( mm_is_hermitian(matcode) || mm_is_symmetric(matcode) )
      ){
    if((ret_code=mm_read_mtx_crd_size(f, &M, &N, &NNZ)) != 0){
      fprintf(stderr, "Error: Could not process size\n");
      exit(1);
    }

    matrix = (double complex *)calloc(NNZ, sizeof(double complex));
    int tmp_row, tmp_col;
    double tmp_real, tmp_imag=0.0;
    if(mm_is_complex(matcode)){
      for(int i=0; i<NNZ; i++){
        fscanf(f, "%d %d %lf %lf\n", &tmp_row, &tmp_col, &tmp_real, &tmp_imag);
        printf("(%d %d %d)\n", tmp_row, tmp_col, tmp_row + tmp_col*(tmp_col-1)/2 - 1);
        matrix[(tmp_row + tmp_col*(tmp_col-1)/2)-1] = tmp_real + tmp_imag*Im_I;
        //tmp_row--;
        //tmp_col--;
        //matrix[tmp_row + tmp_col*N] = tmp_real + tmp_imag*Im_I;
        //matrix[tmp_col + tmp_row*M] = conj(tmp_real + tmp_imag*Im_I);
      }
    }
    else if(mm_is_real(matcode)){
      for(int i=0; i<NNZ; i++){
        fscanf(f, "%d %d %lf\n", &tmp_row, &tmp_col, &tmp_real);
        tmp_row--;
        tmp_col--;
        matrix[tmp_row + tmp_col*(tmp_col+1)/2] = tmp_real;
        matrix[tmp_col + tmp_row*M] = tmp_real;
      }
    }
    else{
      fprintf(stderr, "Error: %s\n", mm_typecode_to_str(matcode));
      fprintf(stderr, "This app require \"real\" or \"complex\"");
      exit(1);
    }
  }
  else{
    fprintf(stderr, "Error: %s\n", mm_typecode_to_str(matcode));
    fprintf(stderr, "This app require \"matrix cordinate complex hermitian\"\n");
    fprintf(stderr, "or\n");
    fprintf(stderr, "This app require \"matrix cordinate real symmetric\"\n");
  }

  char *trans="N"; int m=M,n=N; double complex alpha=1;int ldA=M; int incx=1; double complex beta=0.0; int incy=1;
  double complex *x, *y;
  x = (double complex *)calloc(N, sizeof(double complex));
  y = (double complex *)calloc(N, sizeof(double complex));
  for(int i=0; i<n; i++) x[i]=1.0;
  //zgemv_(trans, &m, &n, &alpha, matrix, &ldA, x, &incx, &beta, y, &incy);
  char *uplo = "U";
  zhpmv_(uplo, &n, &alpha, matrix, x, &incx, &beta, y, &incy);
  for(int i=0; i<n; i++) fprintf(stdout, "%lf %lf\n", creal(y[i]), cimag(y[i]));

  return 0;
}
