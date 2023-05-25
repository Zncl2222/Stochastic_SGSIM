// Copyright 2022 Zncl2222

# include <math.h>

# include "../include/cov_model.h"


void cov_model_init(
    cov_model_t* _cov_model, int bw_l,
    int bw_s, double k_range, double sill, double nugget) {
    _cov_model->bw_l = bw_l;
    _cov_model->bw_s = bw_s;
    _cov_model->bw = bw_l / bw_s;
    _cov_model->k_range = k_range;
    _cov_model->sill = sill;
    _cov_model->nugget = nugget;
}


void cov_model(const double* x, double* cov, int n_dim, cov_model_t* _cov_model) {
    double partial_sill = _cov_model->sill - _cov_model->nugget;
    for (int i = 0; i < n_dim; i++) {
        double factor =
            (1 - exp(-3 * (x[i] * x[i]) / (_cov_model->k_range * _cov_model->k_range)));
        cov[i] = _cov_model->sill - (partial_sill * factor) + _cov_model->nugget;
    }
}

void cov_model2d(const double** x , double* cov, int n_dim, cov_model_t* _cov_model) {
    double partial_sill = _cov_model->sill - _cov_model->nugget;

    for (int i = 0; i < n_dim; i++) {
        for (int j = 0; j < n_dim; j++) {
            double factor =
                (1 - exp(-3 * (x[i][j] * x[i][j]) / (_cov_model->k_range * _cov_model->k_range)));
            cov[n_dim * i + j] = _cov_model->sill - (partial_sill * factor) + _cov_model->nugget;
        }
    }
}
