/**
 * @file sgsim.c
 * @brief Implementation of Sequential Gaussian Simulation (SGSIM) functions.
 *
 * This file contains the implementation of functions for conducting
 * Sequential Gaussian Simulation (SGSIM). It utilizes various libraries and
 * tools for simulations and random number generation.
 *
 * Copyright (c) 2022 Zncl2222
 * License: MIT
 */

# include <stdio.h>
# include <malloc.h>
# include <stdlib.h>
# include <math.h>
# include <float.h>

# include "../include/sgsim.h"
# include "../include/kriging.h"
# include "../include/cov_model.h"
# include "../include/random_tools.h"
# include "../include/matrix_tools.h"
# include "../include/variogram.h"
# include "../include/sort_tools.h"
# include "../lib/c_array.h"

static c_array_int x_grid;
static c_array_double u_array;
static c_array_double sampled;
static c_array_double variogram_array;
static c_array_double sgsim_array;
static sampling_state _sampling;

static int flag;
static int count;
static double epsilon = 1e-6;

void set_sgsim_default(sgsim_t* sgsim, cov_model_t* cov_model) {
    set_cov_model_default(cov_model);
    double boundary_val = pow(cov_model->sill, 0.5) * 4;

    sgsim->z_min = fabs(sgsim->z_min) < epsilon ? -(boundary_val) : sgsim->z_min;
    sgsim->z_max = fabs(sgsim->z_max) < epsilon ? boundary_val : sgsim->z_max;
    sgsim->iteration_limit = sgsim->iteration_limit == 0 ? 10 : sgsim->iteration_limit;

    if (sgsim->if_alloc_memory == 1) {
        unsigned long size = (long)sgsim->x_len * (long)sgsim->realization_numbers;
        sgsim->array = calloc(size, sizeof(double));
    }
}

void sgsim_run(sgsim_t* sgsim, cov_model_t* cov_model, int vario_flag) {
    mt19937_state rng_state;
    mt19937_init(&rng_state, sgsim->randomseed);

    set_sgsim_default(sgsim, cov_model);
    sampling_state_init(&_sampling, sgsim->x_len);

    c_array_init(&variogram_array, cov_model->bw);
    c_array_init(&sgsim_array, sgsim->x_len);

    kriging_param_setting(
        sgsim->x_len, cov_model);  // Initialize parameters

    x_grid.data = arange(sgsim->x_len);
    count = 0;
    int use_cov_cache = 0;
    int error_times = 0;  // The counts of continuosly error
    while (count < sgsim->realization_numbers) {
        printf("Number = %d\n", count);
        _sampling.currlen = 0;
        _sampling.neighbor = 0;
        flag = 0;

        if (cov_model->use_cov_cache == 0 || count == 0) {
            x_grid.data = randompath(x_grid.data, sgsim->x_len, &rng_state);
        } else if (cov_model->use_cov_cache == 0 && count > 0) {
            use_cov_cache = 1;
        }

        for (int i = 0; i < sgsim->x_len; i++) {
            sampling_state_update(&_sampling, x_grid.data[i], i);
            simple_kriging(
                sgsim_array.data,
                &_sampling, &rng_state,
                sgsim->kriging_method,
                use_cov_cache);
            if ((sgsim_array.data[x_grid.data[i]] >= sgsim->z_max)
                 || (sgsim_array.data[x_grid.data[i]] <= sgsim->z_min)) {
                flag++;
                break;
            }
            sgsim->array[x_grid.data[i]+sgsim->x_len*count] = sgsim_array.data[x_grid.data[i]];

            if (_sampling.neighbor < cov_model->max_neighbor) {
                _sampling.neighbor++;
            }

            _sampling.sampled.data[i] = x_grid.data[i];
            _sampling.currlen++;
            if (isfinite(sgsim_array.data[x_grid.data[i]]) == 0) {
                flag++;
            }
        }

        if (flag == 0) {
            save_1darray(sgsim_array.data, sgsim->x_len, "Realizations",
                        "./Realizations/", sgsim->realization_numbers, count);
            if (vario_flag ==1) {
                variogram(
                    sgsim_array.data, variogram_array.data,
                    sgsim->x_len, cov_model->bw, cov_model->bw_s);
                save_1darray(variogram_array.data, cov_model->bw,
                            "Variogram",
                            "./Realizations/Variogram/",
                            sgsim->realization_numbers, count);
            }
            count++;
            error_times = 0;
        } else {
            error_times++;
            if (error_times == sgsim->iteration_limit) {
                fprintf(stderr, "Maximum error occurred. Exiting the program...\n");
                return;
            }
        }
    }
    kriging_memory_free();
    sgsim_memory_free();
}

void sgsim_t_free(sgsim_t* sgsim) {
    free(sgsim->array);
}

static void sgsim_memory_free() {
    c_array_free(&_sampling.sampled);
    c_array_free(&_sampling.u_array);
    c_array_free(&sgsim_array);
    c_array_free(&x_grid);
    c_array_free(&variogram_array);
}
