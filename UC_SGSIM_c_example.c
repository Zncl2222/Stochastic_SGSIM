// Copyright 2022 Zncl2222

# include <stdio.h>
# include <stdlib.h>

# include "./uc_sgsim/c_core/include/sgsim.h"
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

    int mlen = 150;
    int nR = 5;
    int hs = 35;
    int bw = 1;
    int vario_flag = 0;
    int randomseed = 12345;

    sgsim(mlen, nR, hs, bw, 17.32, 1, randomseed, vario_flag);
    PAUSE
    return 0;
}
