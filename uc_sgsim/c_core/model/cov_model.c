// Copyright 2022 Zncl2222

# include <math.h>

# include "../header/cov_model.h"


void Cov_model(double *x, double* cov, int n_dim, double a, double C0) {
    for (int i = 0; i < n_dim; i++) {
        cov[i] = C0 - (C0 * (1 - exp(-3 * pow(x[i], 2) / pow(a, 2))));
    }
}

void Cov_model2d(double **x , double* cov, int n_dim, double a, double C0) {
    for (int i = 0; i < n_dim; i++) {
        for (int j = 0; j < n_dim; j++) {
            cov[n_dim * i + j] =
            C0 - (C0 * (1 - exp(-3 * pow(x[i][j], 2) / pow(a, 2))));
        }
    }
}
