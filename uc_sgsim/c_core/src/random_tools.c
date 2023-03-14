// Copyright 2022 Zncl2222

# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "../include/random_tools.h"
# include "../lib/c_array.h"

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
