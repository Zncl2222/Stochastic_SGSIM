// Copyright 2022 Zncl2222


# include <stdio.h>
# include "../include/cov_model.h"
# include "../include/sgsim.h"
# include "../include/krige.h"
# include "utest.h"

UTEST(test, sgsim) {
    struct sgsim_t sgsim_example;
    sgsim_init(&sgsim_example, 150, 5, 12345, 1);

    struct cov_model_t cov_example;
    cov_model_init(&cov_example, 1, 35, 17.32, 1);

    sgsim_run(&sgsim_example, &cov_example, 0);
    sgsim_run(&sgsim_example, &cov_example, 1);
    sgsim_t_free(&sgsim_example);
}

UTEST_MAIN();
