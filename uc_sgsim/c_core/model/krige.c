// Copyright 2022 Zncl2222

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "..\header\krige.h"
# include "..\header\random_tools.h"
# include "..\header\cov_model.h"
# include "..\header\matrix_tools.h"
# include "..\header\sort_tools.h"

static double range;
static double sill;
static double** array2d_temp;
static double* dist_temp;
static double* distcov_temp;
static double* data_temp;
static double* distcov_temp2;
static double* flatten_temp;
static double* weights;
static double** pdist_temp;
static double** datacov_temp;
static double zvalue, svar;
static float fix;

void Krige_paramsetting(int mlen, double a, double C0) {
    range = a;
    sill = C0;
    dist_temp = (double*)malloc(6 * sizeof(double));
    distcov_temp = (double*)malloc(6 * sizeof(double));
    distcov_temp2 = (double*)malloc(6 * sizeof(double));
    flatten_temp = (double*)malloc(6 * 6 * sizeof(double));
    weights = (double*)malloc(6 * sizeof(double));
    pdist_temp = (double**)malloc(6 * sizeof(double*));
    datacov_temp = (double**)malloc(6 * sizeof(double*));
    array2d_temp = (double**)malloc(mlen * sizeof(double*));

    for (int i = 0; i < mlen; i++) {
        array2d_temp[i] = (double*)malloc(3 * sizeof(double));
    }

    for (int i = 0; i < 6; i++) {
        pdist_temp[i] = (double*)malloc(6 * sizeof(double));
        datacov_temp[i] = (double*)malloc(6 * sizeof(double));
    }
}

void SimpleKrige(double* array, double* sampled, double* u_array, int currlen,
                double unsampled_point, int idx, int neighbor) {
    if (neighbor == 0) {
        array[(int)unsampled_point] = random_normal();
        sampled[idx] = unsampled_point;
    } else {
        int close = 0;

        for (int j = 0; j < currlen; j++) {
            u_array[j] = fabs(sampled[j] - unsampled_point);
            if (u_array[j] < range * 1.732) {
                close++;
            }
        }

        if (close == 0) {
            array[(int)unsampled_point] = random_normal();
            sampled[idx] = unsampled_point;
        } else {
            for (int j = 0; j < currlen; j++) {
                array2d_temp[j][0] = sampled[j];
                array2d_temp[j][1] = array[(int)sampled[j]];
                array2d_temp[j][2] = u_array[j];
            }

            if (neighbor >= 2) {
                // array2d_temp = sort2d(array2d_temp, currlen);
                quicksort2d(array2d_temp, 0, currlen - 1);
            }

            for (int j = 0; j < neighbor; j++) {
                dist_temp[j] = array2d_temp[j][0];
                distcov_temp2[j] = array2d_temp[j][2];
            }

            pdist(dist_temp, pdist_temp, neighbor);
            Cov_model2d(pdist_temp, flatten_temp, neighbor, range, sill);
            matrixform(flatten_temp, datacov_temp, neighbor);
            Cov_model(distcov_temp2, distcov_temp, neighbor, range, sill);

            if (neighbor >= 1)
                LUdecomposition(datacov_temp, distcov_temp, weights, neighbor);

            zvalue = 0;
            svar = 0;
            fix = 0;

            for (int j = 0; j < neighbor; j++) {
                zvalue = zvalue + array2d_temp[j][1] * weights[j];
                svar = svar + distcov_temp[j] * weights[j];
            }

            svar = 1 - svar;
            if (svar < 0)
                svar = 0;
            fix = random_normal() * pow(svar, 0.5);

            array[(int)unsampled_point] = zvalue + fix;
            // printf("VALUE = %lf \n", zvalue + fix);
            // Print_Log1(dist_temp,*array2d_temp,distcov_temp2,currlen,neighbor,unsampled_point);
            // Print_Log2(pdist_temp,datacov_temp,distcov_temp,weights,zvalue,fix,neighbor);
        }
    }
}

void Print_Log1(double* a, double* b, double* c,
                int curr, int n_dim, double u) {
    printf("\nU=%d, current=%d\n\n", (int)u, curr);
    printf("Location = ");
    for (int j = 0; j < n_dim; j++) {
        printf("%lf   ", a[j]);
    }

    printf("\nData =    ");
    for (int j = 0; j < n_dim; j++) {
        printf("%lf   ", b[j]);
    }

    printf("\nDistdiff = ");
    for (int j = 0; j < n_dim; j++) {
        printf("%lf   ", c[j]);
    }
    printf("\n");
}

void Print_Log2(double** a, double** b, double* c, double* d,
                double z_temp, double fix_temp, int n_dim) {
    printf("\nPdist = \n");
    for (int j = 0; j < n_dim; j++) {
        for (int k = 0; k < n_dim; k++) {
            printf("%lf ", a[j][k]);
        }
        printf("\n");
    }
    printf("\n");

    printf("DataCov = \n");
    for (int j = 0; j < n_dim; j++) {
        for (int k = 0; k < n_dim; k++) {
            printf("%lf ", b[j][k]);
        }
        printf("\n");
    }
    printf("\n");

    printf("Distcov = \n");
    for (int j = 0; j < n_dim; j++) {
        printf("%lf   ", c[j]);
    }
    printf("\n");

    printf("Weights = \n");
    for (int j = 0; j < n_dim; j++) {
        printf("%f ", d[j]);
    }

    printf("\n\nFix= %f ", fix_temp);
    printf("\nEstimate=%lf\n", z_temp);
}
