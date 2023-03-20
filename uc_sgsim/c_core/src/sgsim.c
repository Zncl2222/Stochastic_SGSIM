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

static int currentlen;
static int neighbor;
static int flag;
static int count;


void sgsim_init(
    struct sgsim_t* _sgsim,
    int x_grid,
    int realization_numbers,
    int randomseed,
    int if_alloc_memory
) {
    _sgsim->x_grid = x_grid;
    _sgsim->realization_numbers = realization_numbers;
    _sgsim->randomseed = randomseed;
    if (if_alloc_memory == 1) {
        _sgsim->array = malloc(x_grid * realization_numbers * sizeof(double));
    }
}

void sgsim_run(struct sgsim_t* _sgsim, const struct cov_model_t* _cov_model, int vario_flag) {
    mt19937_state rng_state;
    mt19937_init(&rng_state, _sgsim->randomseed);

    c_array_init(&u_array, _sgsim->x_grid);
    c_array_init(&sampled, _sgsim->x_grid);
    c_array_init(&variogram_array, _cov_model->hs);
    c_array_init(&sgsim_array, _sgsim->x_grid);

    Krige_paramsetting(_cov_model->range, _cov_model->sill);  // Initialize parameters

    x_grid.data = arange(_sgsim->x_grid);
    count = 0;
    while (count < _sgsim->realization_numbers) {
        printf("Number = %d\n", count);
        currentlen = 0;
        neighbor = 0;
        flag = 0;

        x_grid.data = randompath(x_grid.data, _sgsim->x_grid, &rng_state);

        for (int i = 0; i < _sgsim->x_grid; i ++) {
            sgsim_array.data[i] = 0;
            sampled.data[i] = 0;
            u_array.data[i] = -1;
        }

        for (int i = 0; i < _sgsim->x_grid; i++) {
            SimpleKrige(sgsim_array.data, sampled.data, u_array.data,
                        currentlen, x_grid.data[i], i, neighbor, &rng_state);
            _sgsim->array[x_grid.data[i]+_sgsim->x_grid*count] = sgsim_array.data[x_grid.data[i]];

            if (neighbor < 8) {
                neighbor++;
            }

            sampled.data[i] = x_grid.data[i];
            currentlen++;
            if (isfinite(sgsim_array.data[x_grid.data[i]]) == 0) {
                flag++;
            }
        }
        count++;

        if (vario_flag == 1)
            variogram(
                sgsim_array.data, variogram_array.data,
                _sgsim->x_grid, _cov_model->hs, _cov_model->bw);

        if (flag > 0) {
            count--;
        } else {
            save_1darray(sgsim_array.data, _sgsim->x_grid, "Realizations",
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

void sgsim_t_free(struct sgsim_t* _sgsim) {
    free(_sgsim->array);
}

static void sgsim_memory_free() {
    c_array_free(&sampled);
    c_array_free(&u_array);
    c_array_free(&sgsim_array);
    c_array_free(&x_grid);
    c_array_free(&variogram_array);
}
