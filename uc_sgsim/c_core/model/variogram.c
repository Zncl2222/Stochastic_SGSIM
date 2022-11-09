// Copyright 2022 Zncl2222

# include <stdio.h>
# include <math.h>
# include <malloc.h>
# include "..\header\variogram.h"
# include "..\header\matrix_tools.h"

void variogram(double* array, double* v , int mlen, int hs, int steps) {
    double Z_temp;
    double* temp;
    double count;
    double** pdist_temp = (double**)malloc(mlen*sizeof(double));

    for (int i = 0; i < mlen; i++) {
        pdist_temp[i] = (double*)malloc(mlen * sizeof(double));
    }

    temp = d_arange(mlen);

    pdist(temp, pdist_temp, mlen);

    for (int i = 0; i < hs; i += steps) {
        Z_temp = 0;
        count = 0;
        for (int j = 0; j < mlen; j++) {
            for (int k = j + 1; k < mlen; k++) {
                if (pdist_temp[j][k] >= i - steps &&
                    pdist_temp[j][k] <= i + steps) {
                    Z_temp = Z_temp+pow((array[j] - array[k]), 2);
                    count += 1;
                }
            }
        }

        if (Z_temp >= 1e-6) {
            v[i] = Z_temp / (2 * count);
        }
    }
}

double variance(double* array, int mlen) {
    double mean = 0;
    double var = 0;
    for (int i = 0; i < mlen; i++) {
        mean = mean + array[i];
    }

    mean = mean / mlen;

    for (int i = 0; i < mlen; i++) {
        var = var +(pow(array[i] - mean, 2));
    }
    var = var / mlen;

    return var;
}
