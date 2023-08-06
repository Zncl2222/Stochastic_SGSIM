// Copyright 2022 Zncl2222

# include <stdio.h>
# include <stdlib.h>

# include "./uc_sgsim/c_core/include/sgsim.h"
# include "./uc_sgsim/c_core/include/cov_model.h"
# if defined(__linux__) || defined(__unix__)
# define PAUSE printf("Press Enter key to continue..."); fgetc(stdin);//NOLINT
# elif _WIN32
# define PAUSE system("PAUSE");
# endif

int main() {
    // you can also set z_min and z_max at sgsim_t. Default value is depends on
    // sill value in cov_model_t
    sgsim_t sgsim_example = {
        .x_len = 150,
        .realization_numbers = 5,
        .randomseed = 12345,
        .kriging_method = 1,
        .if_alloc_memory = 1,
    };

    // you can also set max_negibor at cov_model_t. Defualt value is 4.
    cov_model_t cov_example = {
        .bw_l = 35,
        .bw_s = 1,
        .k_range = 17.32,
        .sill = 1,
        .nugget = 0,
    };

    sgsim_run(&sgsim_example, &cov_example, 0);
    sgsim_t_free(&sgsim_example);
    PAUSE
    return 0;
}
