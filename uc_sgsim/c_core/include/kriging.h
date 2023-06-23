// Copyright 2022 Zncl2222

#ifndef UC_SGSIM_C_CORE_INCLUDE_KRIGING_H_
#define UC_SGSIM_C_CORE_INCLUDE_KRIGING_H_

# include "../lib/c_array.h"
# include "../include/cov_model.h"

typedef struct {
    int neighbor;
    int currlen;
    int idx;
    double unsampled_point;
    c_array_double sampled;
    c_array_double u_array;
} sampling_state;

void sampling_state_init(sampling_state* _sampling, int x_grid_len);

void sampling_state_update(sampling_state* _sampling, double unsampled_point, int idx);

void kriging_param_setting(int x_len, const cov_model_t* _cov_model);

void simple_kriging(double* array, sampling_state* _sampling,
                    mt19937_state* rng_state, int kriging_method);

int find_neighbor(double* array, sampling_state* _sampling, mt19937_state* rng_state);

void matrix_agumented(double** mat, int neighbor);

void kriging_memory_free();

#endif  // UC_SGSIM_C_CORE_INCLUDE_KRIGING_H_
