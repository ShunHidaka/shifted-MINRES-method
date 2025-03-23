#include "function_util.h"


void set_shifts(int *M, double complex **shifts){
  *M = 10;
  (*shifts) = (double complex *)calloc(*M, sizeof(double complex));
  for(int m=0; m<*M; m++){
    (*shifts)[m] = 0.01 * cexp(2 * M_PI * Im_I * (m + 0.5) / *M);
  }
}

void set_rhs(const int N, double complex **rhs){
  (*rhs) = (double complex *)calloc(N, sizeof(double complex));
  for(int i=0; i<N; i++){
    (*rhs)[i] = 1.0;
  }
}

void read_mtx(const char *fname, int *n, double complex **matrix){
  FILE *f;
  int retcode;
  MM_typecode matcode;
  // Open file
  if((f=fopen(fname, "r")) == NULL){
    fprintf(stderr, "Error: Could not found %s\n", fname);
    exit(1);
  }
  // Read Matrix Market banner
  if(mm_read_banner(f, &matcode) != 0){
    fprintf(stderr, "Error: Could not process MM banner.\n");
    exit(1);
  }
  // Read matrix elements
  int M, N, NNZ;
  if( mm_is_matrix(matcode) && mm_is_sparse(matcode) && (mm_is_hermitian(matcode) || mm_is_symmetric(matcode)) ){
    if((retcode=mm_read_mtx_crd_size(f, &M, &N, &NNZ)) != 0){
      fprintf(stderr, "Error: Could not process size\n");
      exit(1);
    }
    *n = N;
    (*matrix) = (double complex *)calloc(N*(N+1)/2, sizeof(double complex));

    int tmp_row, tmp_col;
    double tmp_real, tmp_imag;
    if(mm_is_complex(matcode)){
      while( fscanf(f, "%d %d %lf %lf\n", &tmp_row, &tmp_col, &tmp_real, &tmp_imag) != EOF ){
        (*matrix)[(tmp_row + tmp_col*(tmp_col-1)/2) - 1] = tmp_real + tmp_imag*Im_I;
      }
    }
    else if(mm_is_real(matcode)){
      while( fscanf(f, "%d %d %lf\n", &tmp_row, &tmp_col, &tmp_real) != EOF ){
        (*matrix)[(tmp_row + tmp_col*(tmp_col-1)/2) - 1] = tmp_real;
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
}
