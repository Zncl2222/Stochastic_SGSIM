// Copyright 2022 Zncl2222


# include <stdio.h>
# include "../include/cov_model.h"
# include "../include/sgsim.h"
# include "utest.h"

UTEST(test, sgsim) {
    int mlen = 150;
    int nR = 5;
    int hs = 35;
    int bw = 1;
    int randomseed = 12345;

    sgsim(mlen, nR, hs, bw, 17.32, 1, randomseed, 0);
    sgsim(mlen, nR, hs, bw, 17.32, 1, randomseed, 1);
}

UTEST(test, sgsim_dll) {
    int mlen = 150;
    int nR = 5;
    int hs = 35;
    int bw = 1;
    int randomseed = 1455;
    double* x;
    x = malloc(nR * mlen * sizeof(double));
    sgsim_dll(x, mlen, nR, 17.32, 1, randomseed);
    free(x);
}

UTEST_MAIN();
