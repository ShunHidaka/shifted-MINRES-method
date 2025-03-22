#ifndef FUNCTION_UTIL_H
#define FUNCTION_UTIL_H

// 行列ファイル名
extern char *FNAME;
// メインループの最大値
extern int MAX_ITR;
// BLAS用定数
extern int ZERO;
extern int ONE;
extern double dTMP;
extern double complex cTMP;
// M_PIが未定義の場合定義する
#ifndef M_PI
#define M_PI 3.14159265358979
#endif

// FNAME 引数に応じてを定義する
void set_fname(int c, char *v[]);

// Ax=bのベクトルbを作成する
void set_rhsVector(int N, double complex **b, double *norm);

// M個の純虚数を要素に持つ配列sigmaを作成する
void set_shifts(int *M, double complex **sigma);

// CSR形式による疎行列ベクトル積
void SpMV(const int *A_row, const int *A_col, const double complex *A_ele,
          const double complex *x, double complex *b, int N);


// mtxファイルfnameを開き、コメントアウト部分を読み飛ばし、FILEポインタを返す。
// 行列の行数をrow_sizeに、列数をcol_sizeに、非ゼロ要素数をele_sizeに格納する。
FILE* fopen_mtx(const char *fname, const char *mode,
		int *row_size, int *col_size, int *ele_size);

// CSR形式のデータが格納されたファイル"fname"を開き、
// 行数N, 非ゼロ要素数DATASIZE, CRS形式の行列を返す。
void read_csr(const char *fname,
	      int *N, int *DATASIZE,
	      int **row_ptr, int **col_ind, double complex **element);

#endif
