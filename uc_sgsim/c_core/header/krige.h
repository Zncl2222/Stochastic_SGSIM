// Copyright 2022 Zncl2222
#ifndef UC_SGSIM_C_CORE_HEADER_KRIGE_H_
#define UC_SGSIM_C_CORE_HEADER_KRIGE_H_

void Krige_paramsetting(int X, double a, double C0);

void SimpleKrige(double* array, double* sampled, double* u_array,
                int array_size, double unsampled_point, int idx , int neighbor);

void Print_Log1(double* a, double* b, double* c, int curr, int n_dim, double u);

void Print_Log2(double** a, double** b, double* c, double* d,
                double z_temp, double fix_temp, int n_dim);

#endif  // UC_SGSIM_C_CORE_HEADER_KRIGE_H_
