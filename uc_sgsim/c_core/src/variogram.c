/**
 * @file variogram.c
 * @brief Implementation of Variogram Calculation
 *
 * This source file contains the implementation of variogram calculation functions. Variograms are
 * statistical measures used in geostatistics to quantify the spatial variability of a dataset.
 * The functions in this file calculate the experimental variogram and variance of a given dataset.
 *
 * Copyright (c) 2022-2023 Zncl2222
 * License: MIT
 */

# include <stdio.h>
# include <math.h>
# include <malloc.h>

# include "../include/variogram.h"
# include "../include/matrix_tools.h"
# include "../lib/c_array.h"

void variogram(const double* array, double* v, int mlen, int bw, int bw_s) {
    double z_temp;
    double count;
    c_array_double temp;
    c_matrix_double pdist_temp;
    c_matrix_init(&pdist_temp, mlen, mlen);

    temp.data = d_arange(mlen);

    pdist(temp.data, pdist_temp.data, mlen);

    for (int i = 0; i < bw; i += bw_s) {
        z_temp = 0;
        count = 0;
        double bw_lower = i - bw_s;
        double bw_upper = i + bw_s;

        for (int j = 0; j < mlen; j++) {
            double array_j = array[j];

            for (int k = 0; k < j; k++) {
                double pdist_jk = pdist_temp.data[j][k];

                if (pdist_jk >= bw_lower && pdist_jk <= bw_upper) {
                    double diff = array_j - array[k];
                    z_temp += diff * diff;
                    count += 1;
                }
            }
        }

        if (z_temp >= 1e-6) {
            v[i] = z_temp / (2 * count);
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
        var = var + (pow(array[i] - mean, 2));
    }
    var = var / mlen;

    return var;
}
