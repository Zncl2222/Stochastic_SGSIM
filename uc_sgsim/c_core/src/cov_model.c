/**
 * @file cov_model.c
 * @brief Implementation of Covariance Models
 *
 * This source file contains the implementation of covariance model functions used in geostatistics.
 * These functions provide utilities for setting default covariance model parameters and calculating
 * covariance values based on specified models.
 *
 * Copyright (c) 2022-2023 Zncl2222
 * License: MIT
 */

# include "../include/cov_model.h"

void set_cov_model_default(cov_model_t* _cov_model) {
    _cov_model->max_neighbor = _cov_model->max_neighbor == 0 ? 4 : _cov_model->max_neighbor;
    _cov_model->sill = _cov_model->sill == 0 ? 1 : _cov_model->sill;
    _cov_model->bw = _cov_model->bw_l / _cov_model->bw_s;
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
    int index = 0;
    for (int i = 0; i < n_dim; i++) {
        for (int j = 0; j < n_dim; j++) {
            double factor =
                (1 - exp(-3 * (x[i][j] * x[i][j]) / (_cov_model->k_range * _cov_model->k_range)));
            cov[index] = _cov_model->sill - (partial_sill * factor) + _cov_model->nugget;
            index++;
        }
    }
}
