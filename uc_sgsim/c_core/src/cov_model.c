// Copyright 2022 Zncl2222

# include <math.h>

# include "../include/cov_model.h"


void cov_model_init(
    cov_model_t* _cov_model, int bw,
    int hs, double range, double sill) {
    _cov_model->bw = bw;
    _cov_model->hs = hs;
    _cov_model->range = range;
    _cov_model->sill = sill;
}


void cov_model(const double *x, double* cov, int n_dim, double a, double C0) {
    for (int i = 0; i < n_dim; i++) {
        cov[i] = C0 - (C0 * (1 - exp(-3 * (x[i] * x[i]) / (a * a))));
    }
}

void cov_model2d(const double **x , double* cov, int n_dim, double a, double C0) {
    for (int i = 0; i < n_dim; i++) {
        for (int j = 0; j < n_dim; j++) {
            cov[n_dim * i + j] =
            C0 - (C0 * (1 - exp(-3 * (x[i][j] * x[i][j]) / (a * a))));
        }
    }
}
