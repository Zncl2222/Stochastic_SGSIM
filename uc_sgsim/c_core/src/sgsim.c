// Copyright 2022 Zncl2222


# include <stdio.h>
# include <malloc.h>
# include <stdlib.h>
# include <math.h>

# include "../include/sgsim.h"
# include "../include/krige.h"
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


void sgsim(int X, int nR, int hs, int bw,
        double range, double sill, int randomseed, int vario_flag) {
    mt19937_state rng_state;
    mt19937_init(&rng_state, randomseed);

    c_array_init(&u_array, X);
    c_array_init(&sampled, X);
    c_array_init(&variogram_array, hs);
    c_array_init(&sgsim_array, X);

    Krige_paramsetting(range, sill);  // Initialize parameters

    x_grid.data = arange(X);
    count = 0;
    while (count < nR) {
        printf("Number = %d\n", count);
        currentlen = 0;
        neighbor = 0;
        flag = 0;

        x_grid.data = randompath(x_grid.data, X, &rng_state);

        for (int i = 0; i < X; i ++) {
            sgsim_array.data[i] = 0;
            sampled.data[i] = 0;
            u_array.data[i] = -1;
        }

        for (int i = 0; i < X; i++) {
            SimpleKrige(sgsim_array.data, sampled.data, u_array.data,
                        currentlen, x_grid.data[i], i, neighbor, &rng_state);

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
            variogram(sgsim_array.data, variogram_array.data, X, hs, bw);

        if (flag > 0) {
            count--;
        } else {
            save_1darray(sgsim_array.data, X, "Realizations",
                        "./Realizations/", nR, count);
            if (vario_flag ==1) {
                save_1darray(variogram_array.data, hs,
                            "Variogram", "./Realizations/Variogram/", nR, count);
            }
        }
    }
    krige_memory_free();
    sgsim_memory_free();
    c_array_free(&variogram_array);
}

void sgsim_dll(double* RandomFieldX, int X, int nR,
            double range, double sill, int randomseed) {
    mt19937_state rng_state;
    mt19937_init(&rng_state, randomseed);

    c_array_init(&u_array, X);
    c_array_init(&sampled, X);
    c_array_init(&sgsim_array, X);

    Krige_paramsetting(range, sill);  // Initialize parameters

    x_grid.data = arange(X);
    count = 0;
    while (count < nR) {
        printf("Number = %d\n", count);
        currentlen = 0;
        neighbor = 0;
        flag = 0;

        x_grid.data = randompath(x_grid.data, X, &rng_state);

        for (int i = 0; i < X; i ++) {
            sgsim_array.data[i] = 0;
            sampled.data[i] = 0;
            u_array.data[i] = -1;
        }

        for (int i = 0; i < X; i++) {
            SimpleKrige(sgsim_array.data, sampled.data, u_array.data, currentlen,
                        x_grid.data[i], i, neighbor, &rng_state);
            RandomFieldX[x_grid.data[i]+X*count] = sgsim_array.data[x_grid.data[i]];

            if (neighbor < 8) {
                neighbor++;
            }

            sampled.data[i] = x_grid.data[i];
            currentlen++;
            if (isfinite(sgsim_array.data[x_grid.data[i]]) == 0)
                flag++;
        }
        count++;

        if (flag > 0)
            count--;
    }
    krige_memory_free();
    sgsim_memory_free();
}

void sgsim_memory_free() {
    c_array_free(&sampled);
    c_array_free(&u_array);
    c_array_free(&sgsim_array);
    c_array_free(&x_grid);
}
