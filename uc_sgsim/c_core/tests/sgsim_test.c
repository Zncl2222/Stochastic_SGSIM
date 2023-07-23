// Copyright 2022 Zncl2222


# include <stdio.h>
# include "../include/cov_model.h"
# include "../include/sgsim.h"
# include "../include/kriging.h"
# include "../include/variogram.h"
# include "utest.h"

UTEST(test, sgsim_simple_kriging) {
    sgsim_t sgsim_example;
    sgsim_init(&sgsim_example, 150, 50, 12345, 0, 1);

    cov_model_t cov_example;
    cov_model_init(&cov_example, 35, 1, 17.32, 1, 0);

    sgsim_run(&sgsim_example, &cov_example, 0);
    sgsim_run(&sgsim_example, &cov_example, 1);
    sgsim_t_free(&sgsim_example);
}

UTEST(test, sgsim_simple_kriging2) {
    sgsim_t sgsim_example;
    sgsim_init(&sgsim_example, 300, 16, 987, 0, 1);
    set_sgsim_params(-6, 6, 9);
    cov_model_t cov_example;
    cov_model_init(&cov_example, 35, 1, 25, 1, 0);

    sgsim_run(&sgsim_example, &cov_example, 0);
    sgsim_run(&sgsim_example, &cov_example, 1);
    sgsim_t_free(&sgsim_example);
}

UTEST(test, sgsim_ordinary_kriging) {
    sgsim_t sgsim_example;
    sgsim_init(&sgsim_example, 150, 50, 15, 1, 1);

    cov_model_t cov_example;
    cov_model_init(&cov_example, 35, 1, 17.32, 1, 0);

    sgsim_run(&sgsim_example, &cov_example, 0);
    sgsim_run(&sgsim_example, &cov_example, 1);
    sgsim_t_free(&sgsim_example);
}

UTEST(test, sgsim_ordinary_kriging2) {
    sgsim_t sgsim_example;
    sgsim_init(&sgsim_example, 75, 50, 15, 1, 1);
    set_sgsim_params(-3, 3, 6);
    cov_model_t cov_example;
    cov_model_init(&cov_example, 35, 1, 12, 1, 0);

    sgsim_run(&sgsim_example, &cov_example, 0);
    sgsim_run(&sgsim_example, &cov_example, 1);
    sgsim_t_free(&sgsim_example);
}

UTEST(test, variance) {
    double** arr;
    arr = malloc(20 * sizeof(double));
    for (int i = 0; i < 20; i++) {
        arr[i] = i + i * 2;
    }

    double var = variance(arr, 20);
    free(arr);
}

UTEST_MAIN();
