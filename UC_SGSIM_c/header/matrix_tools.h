# ifndef MATRIX_TOOLS_H
# define MATRIX_TOOLS_H

void LUdecomposition(double** a ,double* b, double* x , int n);

int* arange(int x);
 
double* d_arange(int x);

void pdist(double* x, double** c,int n_dim);

void matrixform(double* x,double**matrix,int n_dim);

void flatten(double** x, int n_dim);

int** matrixReshape(int** mat, int matSize, int* matColSize, int r, int c);

void save_1darray(double* array, int array_size, char* fhead, char* path ,int n_realizations);

# endif