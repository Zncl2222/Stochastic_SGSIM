// Copyright 2022 Zncl2222

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "../include/kriging.h"
# include "../include/random_tools.h"
# include "../include/cov_model.h"
# include "../include/matrix_tools.h"
# include "../include/sort_tools.h"
# include "../lib/c_array.h"

static cov_model_t* model;
static double k_range;
static double estimation;
static double kriging_var;
static double fix;

static c_array_double location;
static c_array_double location_cov;
static c_array_double location_cov2d;
static c_array_double flatten_temp;
static c_array_double weights;
static c_matrix_double distance_mat;
static c_matrix_double pdist_temp;
static c_matrix_double data_cov;

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

void kriging_param_setting(int x_len, const cov_model_t* _cov_model) {
    model = _cov_model;
    k_range = _cov_model->k_range;
    int buffer = _cov_model->max_neighbor + 2;
    c_array_init(&location, buffer);
    c_array_init(&location_cov, buffer);
    c_array_init(&location_cov2d, buffer);
    c_array_init(&weights, buffer);
    c_array_init(&flatten_temp, buffer * buffer);
    c_matrix_init(&pdist_temp, buffer, buffer);
    c_matrix_init(&data_cov, buffer, buffer);
    c_matrix_init(&distance_mat, x_len, 3);
}

void simple_kriging(
    double* array,
    sampling_state* _sampling,
    mt19937_state* rng_state,
    int kriging_method) {
    int has_neighbor = find_neighbor(array, _sampling, rng_state);

    if (has_neighbor == 0) {
        return;
    }

    for (int j = 0; j < _sampling->currlen; j++) {
        distance_mat.data[j][0] = _sampling->sampled.data[j];
        distance_mat.data[j][1] = array[(int)_sampling->sampled.data[j]];
        distance_mat.data[j][2] = _sampling->u_array.data[j];
    }

    if (_sampling->neighbor >= 2) {
        quickselect2d(distance_mat.data, 0, _sampling->currlen - 1, _sampling->neighbor);
    }

    for (int j = 0; j < _sampling->neighbor; j++) {
        location.data[j] = distance_mat.data[j][0];
        location_cov2d.data[j] = distance_mat.data[j][2];
    }

    pdist(location.data, pdist_temp.data, _sampling->neighbor);
    cov_model2d(pdist_temp.data, flatten_temp.data, _sampling->neighbor, model);
    matrixform(flatten_temp.data, data_cov.data, _sampling->neighbor);
    cov_model(location_cov2d.data, location_cov.data, _sampling->neighbor, model);

    if (kriging_method == 1) {
        matrix_agumented(data_cov.data, _sampling->neighbor);
        location_cov.data[_sampling->neighbor] = 1.0;
    }

    int neighbor = kriging_method == 1 ? _sampling->neighbor + 1 : _sampling->neighbor;
    if (_sampling->neighbor >= 1)
        lu_inverse_solver(data_cov.data, location_cov.data, weights.data, neighbor);

    estimation = 0;
    kriging_var = 0;
    fix = 0;

    for (int j = 0; j < _sampling->neighbor; j++) {
        estimation = estimation + distance_mat.data[j][1] * weights.data[j];
        kriging_var = kriging_var + location_cov.data[j] * weights.data[j];
    }

    kriging_var = model->sill - kriging_var;
    if (kriging_var < 0)
        kriging_var = 0;
    fix = mt19937_random_normal(rng_state) * pow(kriging_var, 0.5);

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

void matrix_agumented(double** mat, int neighbor) {
    for (int i = 0; i < neighbor; i++) {
        mat[i][neighbor] = 1;
    }
    for (int i = 0; i <= neighbor; i++) {
        if (i == (neighbor)) {
            mat[neighbor][i] = 0;
        } else {
            mat[neighbor][i] = 1;
        }
    }
}

void kriging_memory_free() {
    c_array_free(&location);
    c_array_free(&location_cov);
    c_array_free(&location_cov2d);
    c_array_free(&flatten_temp);
    c_array_free(&weights);
    c_matrix_free(&distance_mat);
    c_matrix_free(&pdist_temp);
    c_matrix_free(&data_cov);
}
