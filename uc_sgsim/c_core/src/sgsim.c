// Copyright 2022 Zncl2222


# include <stdio.h>
# include <malloc.h>
# include <stdlib.h>
# include <math.h>

# include "../include/sgsim.h"
# include "../include/krige.h"
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

void sgsim_init(
    sgsim_t* _sgsim,
    int x_len,
    int realization_numbers,
    int randomseed,
    int if_alloc_memory
) {
    _sgsim->x_len = x_len;
    _sgsim->realization_numbers = realization_numbers;
    _sgsim->randomseed = randomseed;
    if (if_alloc_memory == 1) {
        unsigned long size = (long)x_len * (long)realization_numbers;
        _sgsim->array = malloc(size * sizeof(double));
    }
}

void sgsim_run(sgsim_t* _sgsim, const cov_model_t* _cov_model, int vario_flag) {
    mt19937_state rng_state;
    mt19937_init(&rng_state, _sgsim->randomseed);

    sampling_state_init(&_sampling, _sgsim->x_len);
    c_array_init(&variogram_array, _cov_model->hs);
    c_array_init(&sgsim_array, _sgsim->x_len);

    krige_param_setting(_cov_model->range, _cov_model->sill);  // Initialize parameters

    x_grid.data = arange(_sgsim->x_len);
    count = 0;
    while (count < _sgsim->realization_numbers) {
        printf("Number = %d\n", count);
        _sampling.currlen = 0;
        _sampling.neighbor = 0;
        flag = 0;

        x_grid.data = randompath(x_grid.data, _sgsim->x_len, &rng_state);

        for (int i = 0; i < _sgsim->x_len; i ++) {
            sgsim_array.data[i] = 0;
            _sampling.sampled.data[i] = 0;
            _sampling.u_array.data[i] = -1;
        }

        for (int i = 0; i < _sgsim->x_len; i++) {
            sampling_state_update(&_sampling, x_grid.data[i], i);
            simple_kriging(sgsim_array.data, &_sampling, &rng_state);
            _sgsim->array[x_grid.data[i]+_sgsim->x_len*count] = sgsim_array.data[x_grid.data[i]];

            if (_sampling.neighbor < 8) {
                _sampling.neighbor++;
            }

            _sampling.sampled.data[i] = x_grid.data[i];
            _sampling.currlen++;
            if (isfinite(sgsim_array.data[x_grid.data[i]]) == 0) {
                flag++;
            }
        }
        count++;

        if (vario_flag == 1)
            variogram(
                sgsim_array.data, variogram_array.data,
                _sgsim->x_len, _cov_model->hs, _cov_model->bw);

        if (flag > 0) {
            count--;
        } else {
            save_1darray(sgsim_array.data, _sgsim->x_len, "Realizations",
                        "./Realizations/", _sgsim->realization_numbers, count);
            if (vario_flag ==1) {
                save_1darray(variogram_array.data, _cov_model->hs,
                            "Variogram",
                            "./Realizations/Variogram/",
                            _sgsim->realization_numbers, count);
            }
        }
    }
    krige_memory_free();
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
