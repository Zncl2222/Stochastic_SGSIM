// Copyright 2022 Zncl2222

# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "../include/random_tools.h"
# include "../include/c_array.h"

# ifndef M_PI
# define M_PI 3.1415926
# endif

int* randompath(int* rpath, int length, mt19937_state* rng_state) {
    int rtemp;

    for (int i = length - 1; i; i--) {
        int rindex = mt19937_generate(rng_state) % i;

        do {
            rtemp = rpath[rindex];
            rpath[rindex] = rpath[i];
            rpath[i] = rtemp;
        }while(0);
    }

    return rpath;
}

float random_normal() {
    float x = rand() / (float)RAND_MAX;
    float y = rand() / (float)RAND_MAX;
    float z = sqrt(-2 * log(x)) * cos(2 * M_PI * y);

    return z;
}
