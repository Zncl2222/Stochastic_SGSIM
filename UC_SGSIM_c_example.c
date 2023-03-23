// Copyright 2022 Zncl2222

# include <stdio.h>
# include <stdlib.h>

# include "./uc_sgsim/c_core/include/sgsim.h"
# include "./uc_sgsim/c_core/include/cov_model.h"
# include "./uc_sgsim/c_core/lib/c_array.h"
# if defined(__linux__) || defined(__unix__)
# define PAUSE printf("Press Enter key to continue..."); fgetc(stdin);//NOLINT
# elif _WIN32
# define PAUSE system("PAUSE");
# endif

int main() {
    /*
        Parameters:
        (int) mlen = x grid of model
        (int) random seed = random seed
        (int) nR = number of realizations
        (int) hs = lag tolal distance (for variogram calculation)
        (int) bw = steps of lag distance
        (int) vario_flag => 0 won't calulate the variogram.
    */

    sgsim_t sgsim_example;
    // Initialize with function
    sgsim_init(&sgsim_example, 150, 5, 12345, 1);
    // Initialize with assignment
    // sgsim.x_grid = 150;
    // sgsim_example.realization_numbers = 5;
    // sgsim_example.randomseed = 12345;
    // sgsim_example.if_alloc_memory = 0;
    // sgsim_example.array = malloc(150 * 5 * sizeof(double));

    cov_model_t cov_example;
    // Initialize with function
    cov_model_init(&cov_example, 1, 35, 17.32, 1);
    // Initialize with assignment
    // cov_example.bw = 1;
    // cov_example.hs = 35;
    // cov_example.range = 17.32;
    // cov_example.sill = 1;

    sgsim_run(&sgsim_example, &cov_example, 0);
    sgsim_t_free(&sgsim_example);
    PAUSE
    return 0;
}
