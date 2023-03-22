// Copyright 2022 Zncl2222

#ifndef UC_SGSIM_C_CORE_INCLUDE_KRIGE_H_
#define UC_SGSIM_C_CORE_INCLUDE_KRIGE_H_

# include "../lib/c_array.h"

struct sampling_state {
    int neighbor;
    int currlen;
    int idx;
    double unsampled_point;
    c_array_double sampled;
    c_array_double u_array;
};

void sampling_state_init(struct sampling_state* _sampling, int x_grid_len);

void sampling_state_update(struct sampling_state* _sampling, double unsampled_point, int idx);

void krige_param_setting(double a, double C0);

void simple_kriging(double* array, struct sampling_state* _sampling, mt19937_state* rng_state);

int find_neighbor(double* array, struct sampling_state* _sampling, mt19937_state* rng_state);

void krige_memory_free();

#endif  // UC_SGSIM_C_CORE_INCLUDE_KRIGE_H_
