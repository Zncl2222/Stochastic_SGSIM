// Copyright 2022 Zncl2222


# include <stdio.h>
# include <malloc.h>
# include <stdlib.h>
# include <math.h>

# include "..\header\sgsim.h"
# include "..\header\krige.h"
# include "..\header\random_tools.h"
# include "..\header\matrix_tools.h"
# include "..\header\variogram.h"
# include "..\header\sort_tools.h"

int* x_grid;
int currentlen;
int neighbor;
int flag;

double* sampled;
double* u_array;
double* variogram_array;
double* sgsim_array;

int count;

void sgsim(int X, int nR, int hs, int bw,
        double range, double sill, int vario_flag) {
    u_array = (double*)malloc(X * sizeof(double));
    x_grid = (int*)malloc(X * sizeof(int));
    sampled = (double*)malloc(X * sizeof(double));

    variogram_array = (double*)malloc(hs * sizeof(double));
    sgsim_array = (double*)malloc(X * sizeof(double));

    Krige_paramsetting(X, range, sill);  // Initialize parameters

    x_grid = arange(X);
    count = 0;
    while (count < nR) {
        printf("Number = %d\n", count);
        currentlen = 0;
        neighbor = 0;
        flag = 0;

        x_grid = randompath(x_grid, X);

        for (int i = 0; i < X; i ++) {
            sgsim_array[i] = 0;
            sampled[i] = 0;
            u_array[i] = -1;
        }

        for (int i = 0; i < X; i++) {
            SimpleKrige(sgsim_array, sampled, u_array,
                        currentlen, x_grid[i], i, neighbor);

            if (neighbor < 5) {
                neighbor++;
            }

            sampled[i] = x_grid[i];
            currentlen++;
            if (isfinite(sgsim_array[(int)x_grid[i]]) == 0) {
                flag++;
            }
        }
        count++;

        if (vario_flag == 1)
            variogram(sgsim_array, variogram_array, X, hs, bw);

        if (flag > 0) {
            count--;
        } else {
            save_1darray(sgsim_array, X, "Realizations",
                        "./Realizations/" , count);
            if (vario_flag ==1) {
                save_1darray(variogram_array, hs,
                            "Variogram", "./Realizations/Variogram/", count);
            }
        }
    }
}

void sgsim_dll(double* RandomFieldX, int X, int nR,
            double range, double sill) {
    double* sgsim_array = (double*)malloc(X * sizeof(double));
    u_array = (double*)malloc(X * sizeof(double));
    x_grid = (int*)malloc(X * sizeof(int));
    sampled = (double*)malloc(X * sizeof(double));

    count = 0;
    Krige_paramsetting(X, range, sill);  // Initialize parameters
    x_grid = arange(X);

    while (count < nR) {
        printf("NUmber = %d\n", count);
        currentlen = 0;
        neighbor = 0;
        flag = 0;

        x_grid = randompath(x_grid, X);

        for (int i = 0; i < X; i ++) {
            sgsim_array[i] = 0;
            sampled[i] = 0;
            u_array[i] = -1;
        }

        for (int i = 0; i < X; i++) {
            SimpleKrige(sgsim_array, sampled, u_array, currentlen,
                        x_grid[i], i, neighbor);
            RandomFieldX[(int)x_grid[i]+X*count] = sgsim_array[(int)x_grid[i]];

            if (neighbor < 5)
                neighbor++;

            sampled[i] = x_grid[i];
            currentlen++;
            if (isfinite(sgsim_array[(int)x_grid[i]]) == 0)
                flag++;
        }
        count++;

        if (flag > 0)
            count--;
    }
}
