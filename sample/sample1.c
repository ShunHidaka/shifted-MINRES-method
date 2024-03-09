#include "../c/sminres_solver.h"
char *FNAME = "/home/stern/multi/data/PPE3594.csr";
//char *FNAME = "/home/stern/multi/data/VCNT900h.csr";
int M = 10;
FILE *fopen_mtx(const char *fname, const char *mode, int *row_size, int *col_size, int *ele_size);
void read_csr(const char *fname, int *N, int *DATASIZE, int **row_ptr, int **col_ind, double complex **element);

int main(void){
  int i;

  /* read R-Symmetric of C-Hermitian Matrix 'A' */
  int N, DATASIZE;
  int *A_row, *A_col;
  double complex *A_ele;
  read_csr(FNAME, &N, &DATASIZE, &A_row, &A_col, &A_ele);

  /* generate shift points */
  double complex *sigma;
  sigma = (double complex *)calloc(M, sizeof(double complex));
  for(i=0; i<M; i++)
    sigma[i] = 0.01 * cexp( 2 * M_PI * I * (i + 0.5) / M );

  /* generate rhs-vector */
  double complex * b;
  b = (double complex *)calloc(N, sizeof(double complex));
  for(i=0; i<N; i++)
    b[i] = 1.0;

  int result;
  int itermax = 100000;
  double threshold = 1e-13;
  double complex *x;
  int status[3];
  x = (double complex *)calloc(M*N, sizeof(double complex));

  fprintf(stdout, "%s\n", FNAME);
  status[0]=1; status[1]=0; status[2]=20;
  result = solve_sminres(sigma, M, A_row,A_col,A_ele, x, b, N, itermax, threshold, status);

  if(result == 0)
    printf("fail\n");
  else
    printf("success\n");
  
  return 0;
}

FILE *fopen_mtx(const char *fname, const char *mode, int *row_size, int *col_size, int *ele_size)
{
  FILE *fp = NULL;
  fp = fopen(fname, mode);
  if(fp == NULL){
    fprintf(stderr, "Can not open file : %s\n", fname);
    exit(1);
  }
  char chr;
  while( (chr = fgetc(fp)) != EOF && (chr=='%' || chr=='#') )
    while( (chr = fgetc(fp)) != EOF )
      if(chr == '\n') break;
  fseek(fp, -sizeof(char), SEEK_CUR);
  int num, tmp1, tmp2, tmp3;
  num = fscanf(fp, "%d %d %d", &tmp1, &tmp2, &tmp3);
  if(num != 3){ fprintf(stderr, "fopen err\n"); exit(1);}
  if(col_size != NULL) *row_size = tmp1;
  if(row_size != NULL) *col_size = tmp2;
  if(ele_size != NULL) *ele_size = tmp3;
  return fp;
}
void read_csr(const char *fname, int *N, int *DATASIZE,
              int **row_ptr, int **col_ind, double complex **element)
{
  FILE *fp;
  int rsize, csize, esize;
  fp = fopen_mtx(fname, "r", &rsize, &csize, &esize);
  *N = rsize-1; *DATASIZE = esize;

  *row_ptr = (int *)calloc(rsize, sizeof(int));
  *col_ind = (int *)calloc(csize, sizeof(int));
  *element = (double complex *)calloc(esize, sizeof(double complex));

  int row, col, i=0;
  double real, imag=0;
  while( fscanf(fp, "%d %d %lf %lf", &row, &col, &real, &imag) != EOF ){
    if(i < rsize) (*row_ptr)[i] = row;
    if(i < csize) (*col_ind)[i] = col;
    if(i < esize) (*element)[i] = CMPLX(real, imag);//real + imag*I;
    i = i + 1;
  }
  fclose(fp);
}

/*
gcc sample1.c ../c/sminres_solver.c -lm -lgfortran -lblas -llapack
*/
