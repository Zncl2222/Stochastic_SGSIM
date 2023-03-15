// Copyright 2022 Zncl2222

# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <string.h>
# include <errno.h>
# include "../include/matrix_tools.h"
# ifdef __WIN32__
# include <io.h>
# elif defined(__linux__) || defined(__unix__)
# include <fcntl.h>
# include <sys/io.h>
# include <sys/stat.h>
struct stat st = {0};
# endif

void LUdecomposition(double** a , double* b, double* x , int n) {
    int i = 0, j = 0, k = 0;
    double l[10][10], u[10][10];
    double c[10];

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (j < i) {
                l[j][i] = 0;
            } else {
                l[j][i] = a[j][i];
                for (k = 0; k < i; k++) {
                    l[j][i] = l[j][i] - l[j][k] * u[k][i];
                }
            }
        }
        for (j = 0; j < n; j++) {
            if (j < i) {
                u[i][j] = 0;
            } else if (j == i) {
                u[i][j] = 1;
            } else {
                u[i][j] = a[i][j] / l[i][i];
                for (k = 0; k < i; k++) {
                    u[i][j] = u[i][j] - ((l[i][k] * u[k][j]) / l[i][i]);
                }
            }
        }
    }
    // Solve L(Ux)=b, assume Ux=c
    c[0] = b[0] / l[0][0];

    for (i = 1; i < n; i++) {
        c[i] = b[i];

        for (j = 0; j < i; j++) {
            c[i] = c[i] - l[i][j] * c[j];
        }

        c[i] = c[i] / l[i][i];
    }

    x[n] = c[n];

    for (i = n - 1; i >= 0; i--) {
        x[i] = c[i];

        for (j = i + 1; j < n; j++) {
            x[i] -= u[i][j] * x[j];
        }
    }
}

int* arange(int x) {
    int* res;
    res = (int*)malloc(x*sizeof(int));

    for (int i = 0; i < x; i++) {
        res[i] = i;
    }
    return res;
}

double* d_arange(int x) {
    double* space;
    space = (double*)malloc(x*sizeof(double));

    for (int i = 0; i < x; i++) {
        space[i] = i;
    }
    return space;
}

void pdist(const double* x, double** c, int n_dim) {
    for (int i = 0; i < n_dim; i++) {
        for (int j = 0; j < n_dim; j++) {
            c[i][j] = fabs(x[j] - x[i]);
        }
    }
}

void matrixform(const double* x, double**matrix, int n_dim) {
    for (int i = 0; i < n_dim; i++) {
        for (int j = 0; j < n_dim; j++) {
            matrix[i][j] = x[n_dim*i+j];
        }
    }
}

void save_1darray(double* array, int array_size,
                char* fhead, char* path , int total_n, int curr_n) {
    FILE *output;
    int num_digits = (int)ceil(log10(total_n)) + 1;
    char format[200];

    snprintf(format, sizeof(format), "%s%s%%0%dd.txt", path, fhead, num_digits);

    #ifdef __linux__
    if (mkdir(path, 0777 | O_EXCL) == -1) {
        if (errno != EEXIST) {
            printf("Folder already exist!");
        }
    }
    #elif defined(__WIN32__)
    if (_mkdir(path) == -1) {
        if (errno != EEXIST) {
            printf("Folder already exist!");
        }
    }
    #endif

    char filename[200];
    snprintf(filename, sizeof(filename), format, curr_n);

    output = fopen(filename, "w");
    if (output == NULL) {
        perror("Failed to open the file");
        exit(1);
    }

    for (int i = 0; i < array_size; i++) {
        fprintf(output, "%d\t%.10f\n", i, array[i]);
    }

    fclose(output);
}
