// Copyright 2022 Zncl2222

# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <string.h>
# include <io.h>

# include "..\header\matrix_tools.h"

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

void pdist(double* x, double** c, int n_dim) {
    for (int i = 0; i < n_dim; i++) {
        for (int j = 0; j < n_dim; j++) {
            c[i][j] = fabs(x[j] - x[i]);
        }
    }
}

void matrixform(double* x, double**matrix, int n_dim) {
    for (int i = 0; i < n_dim; i++) {
        for (int j = 0; j < n_dim; j++) {
            matrix[i][j] = x[n_dim*i+j];
        }
    }
}

void flatten(double** x, int n_dim) {
    double* flat;

    flat = (double*)malloc(n_dim*n_dim*sizeof(double));

    for (int i = 0; i < n_dim; i++) {
        for (int j = 0; j < n_dim; j++) {
            flat[n_dim * i + j] = x[i][j];
        }
    }
}

int** matrixReshape(int** mat, int matSize, int* matColSize, int r, int c) {
    int matCol = matColSize[0];

    if (matSize*matCol != r * c || (matSize == r && matCol == c)) {
        printf("(matrixReshape)ERROR: ");
        printf("Matrix can't reshape to the target shape.\n");
        return mat;
    }
    int** res = (int**)malloc(r*sizeof(int*));

    for (int i = 0; i < r; i++) {
        res[i] = (int*)malloc(c*sizeof(int));
    }

    for (int i = 0; i < r * c; i++) {
        res[i / c][i % c] = mat[i / matCol][i % matCol];
    }

    return res;
}

void save_1darray(double* array, int array_size,
                char* fhead, char* path , int n_realizations) {
    FILE *output;

    char n1[] = "000", n2[] = "00", n3[] = "0";
    char ftail[] =".txt";
    char number1[15];
    char fullname[200];
    int max_len =  strlen(path) + strlen(ftail) + strlen(fhead) + 12 + 4;

    mkdir(path);

    snprintf(number1, sizeof(number1), "%d", n_realizations);

    memset(fullname, '\0', max_len);
    snprintf(fullname, max_len - strlen(path), "%s", path);
    snprintf(fullname + strlen(fullname),
            max_len - strlen(fhead) + 1, "%s", fhead);

    if (n_realizations < 10)
        snprintf(fullname + strlen(fullname),
                max_len - strlen(n1) + 1, "%s", n1);
    else if (n_realizations < 100 && n_realizations >= 10)
        snprintf(fullname + strlen(fullname),
                max_len - strlen(n2) + 1, "%s", n2);
    else if (n_realizations < 1000 && n_realizations >= 100)
        snprintf(fullname + strlen(fullname),
                max_len - strlen(n3) + 1, "%s", n3);

    snprintf(fullname + strlen(fullname),
            max_len - strlen(number1), "%s", number1);
    snprintf(fullname + strlen(fullname), max_len - strlen(ftail), "%s", ftail);

    output = fopen(fullname, "w");

    for (int i = 0; i < array_size; i++) {
        fprintf(output, "%d\t%.10f\n", i, array[i]);
    }

    fclose(output);
}
