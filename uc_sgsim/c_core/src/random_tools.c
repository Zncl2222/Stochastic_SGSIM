/**
 * @file random_tools.c
 * @brief Implementation of Randomization Tools
 *
 * This source file contains the implementation of randomization tools, including functions for
 * shuffling arrays and generating random paths. These functions are essential for various
 * randomization tasks.
 *
 * Copyright (c) 2022-2023 Zncl2222
 * License: MIT
 */

# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "../include/random_tools.h"
# include "../c_array_tools/src/c_array.h"

# ifndef M_PI
# define M_PI 3.1415926
# endif

int* randompath(int* rpath, int length, mt19937_state* rng_state) {
    int rtemp;

    for (int i = length - 1; i; i--) {
        int rindex = mt19937_generate(rng_state) % i;
        rtemp = rpath[rindex];
        rpath[rindex] = rpath[i];
        rpath[i] = rtemp;
    }

    return rpath;
}
