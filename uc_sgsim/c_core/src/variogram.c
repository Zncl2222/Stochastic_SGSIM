// Copyright 2022 Zncl2222

# include <stdio.h>
# include <math.h>
# include <malloc.h>

# include "../include/variogram.h"
# include "../include/matrix_tools.h"
# include "../lib/c_array.h"

void variogram(const double* array, double* v , int mlen, int hs, int steps) {
    double Z_temp;
    double count;
    c_array(double) temp;
    c_matrix(double) pdist_temp;
    c_matrix_init(&pdist_temp, mlen, mlen);

    temp.data = d_arange(mlen);

    pdist(temp.data, pdist_temp.data, mlen);

    for (int i = 0; i < hs; i += steps) {
        Z_temp = 0;
        count = 0;
        for (int j = 0; j < mlen; j++) {
            for (int k = j + 1; k < mlen; k++) {
                if (pdist_temp.data[j][k] >= i - steps &&
                    pdist_temp.data[j][k] <= i + steps) {
                    Z_temp = Z_temp+pow((array[j] - array[k]), 2);
                    count += 1;
                }
            }
        }

        if (Z_temp >= 1e-6) {
            v[i] = Z_temp / (2 * count);
        }
    }
    c_array_free(&temp);
    c_matrix_free(&pdist_temp);
}

double variance(const double* array, int mlen) {
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
