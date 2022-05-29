# ifndef KRIGE_H
# define KRIGE_H

void Krige_paramsetting(int X, double a, double C0);

void SimpleKrige(double* array, double* sampled, double* u_array, int array_size, double unsampled_point, int idx ,int neighbor, int randomseed);

void Print_Log1(double* a,double* b,double* c, int curr, int n_dim, double u);

void Print_Log2(double** a,double** b,double* c,double* d,double z_temp, double fix_temp, int n_dim);

# endif