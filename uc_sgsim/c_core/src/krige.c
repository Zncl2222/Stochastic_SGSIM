// Copyright 2022 Zncl2222

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "../include/krige.h"
# include "../include/random_tools.h"
# include "../include/cov_model.h"
# include "../include/matrix_tools.h"
# include "../include/sort_tools.h"
# include "../lib/c_array.h"

static cov_model_t* model;
static double k_range;
static double estimation;
static double krige_var;
static double fix;

static c_array_double location;
static c_array_double loc_cov;
static c_array_double data_temp;
static c_array_double loc_cov2;
static c_array_double flatten_temp;
static c_array_double weights;
static c_matrix_double array2d_temp;
static c_matrix_double pdist_temp;
static c_matrix_double datacov;

void sampling_state_init(sampling_state* _sampling, int x_grid_len) {
    _sampling->neighbor = 0;
    _sampling->currlen = 0;
    _sampling->idx = 0;
    _sampling->unsampled_point = 0;

    c_array_init(&_sampling->sampled, x_grid_len);
    c_array_init(&_sampling->u_array, x_grid_len);
}

void sampling_state_update(sampling_state* _sampling, double unsampled_point, int idx) {
    _sampling->unsampled_point = unsampled_point;
    _sampling->idx = idx;
}

void krige_param_setting(int x_len, const cov_model_t* _cov_model) {
    model = _cov_model;
    k_range = _cov_model->k_range;
    c_array_init(&location, 10);
    c_array_init(&loc_cov, 10);
    c_array_init(&loc_cov2, 10);
    c_array_init(&weights, 10);
    c_array_init(&flatten_temp, 100);
    c_array_init(&data_temp, 10);
    c_matrix_init(&pdist_temp, 10, 10);
    c_matrix_init(&datacov, 10, 10);
    c_matrix_init(&array2d_temp, x_len, 3);
}

void simple_kriging(double* array, sampling_state* _sampling, mt19937_state* rng_state) {
    int has_neighbor = find_neighbor(array, _sampling, rng_state);

    if (has_neighbor == 0) {
        return;
    }

    for (int j = 0; j < _sampling->currlen; j++) {
        array2d_temp.data[j][0] = _sampling->sampled.data[j];
        array2d_temp.data[j][1] = array[(int)_sampling->sampled.data[j]];
        array2d_temp.data[j][2] = _sampling->u_array.data[j];
    }

    if (_sampling->neighbor >= 2) {
        quicksort2d(array2d_temp.data, 0, _sampling->currlen - 1);
    }

    for (int j = 0; j < _sampling->neighbor; j++) {
        location.data[j] = array2d_temp.data[j][0];
        loc_cov2.data[j] = array2d_temp.data[j][2];
    }

    pdist(location.data, pdist_temp.data, _sampling->neighbor);
    cov_model2d(pdist_temp.data, flatten_temp.data, _sampling->neighbor, model);
    matrixform(flatten_temp.data, datacov.data, _sampling->neighbor);
    cov_model(loc_cov2.data, loc_cov.data, _sampling->neighbor, model);

    if (_sampling->neighbor >= 1)
        lu_inverse_solver(datacov.data, loc_cov.data, weights.data, _sampling->neighbor);

    estimation = 0;
    krige_var = 0;
    fix = 0;

    for (int j = 0; j < _sampling->neighbor; j++) {
        estimation = estimation + array2d_temp.data[j][1] * weights.data[j];
        krige_var = krige_var + loc_cov.data[j] * weights.data[j];
    }

    krige_var = model->sill - krige_var;
    if (krige_var < 0)
        krige_var = 0;
    fix = mt19937_random_normal(rng_state) * pow(krige_var, 0.5);

    array[(int)_sampling->unsampled_point] = estimation + fix;
}

int find_neighbor(double* array, sampling_state* _sampling,
                  mt19937_state* rng_state) {
    if (_sampling->neighbor == 0) {
        array[(int)_sampling->unsampled_point] = mt19937_random_normal(rng_state) * model->sill;
        _sampling->sampled.data[_sampling->idx] = _sampling->unsampled_point;
        return 0;
    }
    int close = 0;

    for (int j = 0; j < _sampling->currlen; j++) {
        _sampling->u_array.data[j] = fabs(_sampling->sampled.data[j] - _sampling->unsampled_point);
        if (_sampling->u_array.data[j] < k_range * 1.732) {
            close++;
        }
    }

    if (close == 0) {
        array[(int)_sampling->unsampled_point] = mt19937_random_normal(rng_state) * model->sill;
        _sampling->sampled.data[_sampling->idx] = _sampling->unsampled_point;
        return 0;
    }

    return 1;
}

void krige_memory_free() {
    c_array_free(&location);
    c_array_free(&loc_cov);
    c_array_free(&data_temp);
    c_array_free(&loc_cov2);
    c_array_free(&flatten_temp);
    c_array_free(&weights);
    c_matrix_free(&array2d_temp);
    c_matrix_free(&pdist_temp);
    c_matrix_free(&datacov);
}
