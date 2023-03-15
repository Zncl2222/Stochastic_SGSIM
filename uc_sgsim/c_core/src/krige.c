// Copyright 2022 Zncl2222

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "../include/krige.h"
# include "../include/random_tools.h"
# include "../include/cov_model.h"
# include "../include/matrix_tools.h"
# include "../include/sort_tools.h"
# include "../lib/c_array.h"

static double range;
static double sill;
static double estimation;
static double krige_var;
static double fix;

static c_array(double) location;
static c_array(double) loc_cov;
static c_array(double) data_temp;
static c_array(double) loc_cov2;
static c_array(double) flatten_temp;
static c_array(double) weights;
static c_matrix(double) array2d_temp;
static c_matrix(double) pdist_temp;
static c_matrix(double) datacov;


void Krige_paramsetting(double a, double C0) {
    range = a;
    sill = C0;
    c_array_init(&location, 10);
    c_array_init(&loc_cov, 10);
    c_array_init(&loc_cov2, 10);
    c_array_init(&weights, 10);
    c_array_init(&flatten_temp, 10 * 10);
    c_array_init(&data_temp, 10);
    c_matrix_init(&pdist_temp, 10, 10);
    c_matrix_init(&datacov, 10, 10);
    c_matrix_init(&array2d_temp, 150, 3);
}

void SimpleKrige(double* array, double* sampled, double* u_array, int currlen,
                double unsampled_point, int idx, int neighbor, mt19937_state* rng_state) {
    int has_neighbor = find_neighbor(array, sampled, u_array, currlen,
                                    unsampled_point, idx, neighbor, rng_state);

    if (has_neighbor == 0) {
        return;
    }

    for (int j = 0; j < currlen; j++) {
        array2d_temp.data[j][0] = sampled[j];
        array2d_temp.data[j][1] = array[(int)sampled[j]];
        array2d_temp.data[j][2] = u_array[j];
    }

    if (neighbor >= 2) {
        quicksort2d(array2d_temp.data, 0, currlen - 1);
    }

    for (int j = 0; j < neighbor; j++) {
        location.data[j] = array2d_temp.data[j][0];
        loc_cov2.data[j] = array2d_temp.data[j][2];
    }

    pdist(location.data, pdist_temp.data, neighbor);
    Cov_model2d(pdist_temp.data, flatten_temp.data, neighbor, range, sill);
    matrixform(flatten_temp.data, datacov.data, neighbor);
    Cov_model(loc_cov2.data, loc_cov.data, neighbor, range, sill);

    if (neighbor >= 1)
        LUdecomposition(datacov.data, loc_cov.data, weights.data, neighbor);

    estimation = 0;
    krige_var = 0;
    fix = 0;

    for (int j = 0; j < neighbor; j++) {
        estimation = estimation + array2d_temp.data[j][1] * weights.data[j];
        krige_var = krige_var + loc_cov.data[j] * weights.data[j];
    }

    krige_var = 1 - krige_var;
    if (krige_var < 0)
        krige_var = 0;
    fix = mt19937_random_normal(rng_state) * pow(krige_var, 0.5);

    array[(int)unsampled_point] = estimation + fix;
}

int find_neighbor(double* array, double* sampled, double* u_array, int currlen,
                double unsampled_point, int idx, int neighbor, mt19937_state* rng_state) {
    if (neighbor == 0) {
        array[(int)unsampled_point] = mt19937_random_normal(rng_state);
        sampled[idx] = unsampled_point;
        return 0;
    }
    int close = 0;

    for (int j = 0; j < currlen; j++) {
        u_array[j] = fabs(sampled[j] - unsampled_point);
        if (u_array[j] < range * 1.732) {
            close++;
        }
    }

    if (close == 0) {
        array[(int)unsampled_point] = mt19937_random_normal(rng_state);
        sampled[idx] = unsampled_point;
        return 0;
    }

    return 1;
}

void krige_memory_free() {
    c_array_free(&location);
    c_array_free(&loc_cov);
    c_array_free(&data_temp);
    c_array_free(&loc_cov2);
    c_array_free(&flatten_temp);
    c_array_free(&weights);
    c_matrix_free(&array2d_temp);
    c_matrix_free(&pdist_temp);
    c_matrix_free(&datacov);
}
