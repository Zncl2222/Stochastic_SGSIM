// Copyright 2022 Zncl2222
#ifndef UC_SGSIM_C_CORE_INCLUDE_MATRIX_TOOLS_H_
#define UC_SGSIM_C_CORE_INCLUDE_MATRIX_TOOLS_H_

void lu_decomposition(double** mat, double** l, double** u, int n);

void lu_inverse_solver(double** mat, const double* array, double* result, int n);

int* arange(int x);

double* d_arange(int x);

void pdist(const double* x, double** c, int n_dim);

void matrixform(const double* x, double**matrix, int n_dim);

void save_1darray(const double* array, int array_size,
                char* fhead, char* path, int total_n, int curr_n);

#endif  // UC_SGSIM_C_CORE_INCLUDE_MATRIX_TOOLS_H_
