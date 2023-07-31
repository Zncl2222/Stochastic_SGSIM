// Copyright 2022 Zncl2222


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

void set_sgsim_default(sgsim_t* _sgsim, cov_model_t* _cov_model) {
    set_cov_model_default(_cov_model);
    double boundary_val = pow(_cov_model->sill, 0.5) * 4;
    _sgsim->z_min = _sgsim->z_min == 0 ? -(boundary_val) : _sgsim->z_min;
    _sgsim->z_max = _sgsim->z_max == 0 ? (boundary_val) : _sgsim->z_max;
    unsigned long size = (long)_sgsim->x_len * (long)_sgsim->realization_numbers;
    _sgsim->array = calloc(size, sizeof(double));
}

void sgsim_run(sgsim_t* _sgsim, cov_model_t* _cov_model, int vario_flag) {
    mt19937_state rng_state;
    mt19937_init(&rng_state, _sgsim->randomseed);
    set_sgsim_default(_sgsim, _cov_model);

    sampling_state_init(&_sampling, _sgsim->x_len);
    c_array_init(&variogram_array, _cov_model->bw);
    c_array_init(&sgsim_array, _sgsim->x_len);

    kriging_param_setting(
        _sgsim->x_len, _cov_model);  // Initialize parameters

    x_grid.data = arange(_sgsim->x_len);
    count = 0;
    while (count < _sgsim->realization_numbers) {
        printf("Number = %d\n", count);
        _sampling.currlen = 0;
        _sampling.neighbor = 0;
        flag = 0;

        x_grid.data = randompath(x_grid.data, _sgsim->x_len, &rng_state);
        for (int i = 0; i < _sgsim->x_len; i++) {
            sampling_state_update(&_sampling, x_grid.data[i], i);
            simple_kriging(sgsim_array.data, &_sampling, &rng_state, _sgsim->kriging_method);
            if ((sgsim_array.data[x_grid.data[i]] >= _sgsim->z_max)
                 || (sgsim_array.data[x_grid.data[i]] <= _sgsim->z_min)) {
                flag++;
                break;
            }
            _sgsim->array[x_grid.data[i]+_sgsim->x_len*count] = sgsim_array.data[x_grid.data[i]];

            if (_sampling.neighbor < _cov_model->max_neighbor) {
                _sampling.neighbor++;
            }

            _sampling.sampled.data[i] = x_grid.data[i];
            _sampling.currlen++;
            if (isfinite(sgsim_array.data[x_grid.data[i]]) == 0) {
                flag++;
            }
        }

        if (flag == 0) {
            save_1darray(sgsim_array.data, _sgsim->x_len, "Realizations",
                        "./Realizations/", _sgsim->realization_numbers, count);
            if (vario_flag ==1) {
                variogram(
                    sgsim_array.data, variogram_array.data,
                    _sgsim->x_len, _cov_model->bw, _cov_model->bw_s);
                save_1darray(variogram_array.data, _cov_model->bw,
                            "Variogram",
                            "./Realizations/Variogram/",
                            _sgsim->realization_numbers, count);
            }
            count++;
        }
    }
    kriging_memory_free();
    sgsim_memory_free();
}

void sgsim_t_free(sgsim_t* _sgsim) {
    free(_sgsim->array);
}

static void sgsim_memory_free() {
    c_array_free(&_sampling.sampled);
    c_array_free(&_sampling.u_array);
    c_array_free(&sgsim_array);
    c_array_free(&x_grid);
    c_array_free(&variogram_array);
}
