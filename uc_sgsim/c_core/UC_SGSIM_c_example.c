// Copyright 2022 Zncl2222

# include <stdio.h>
# include <stdlib.h>

# include "header\matrix_tools.h"
# include "header\variogram.h"
# include "header\krige.h"
# include "header\cov_model.h"
# include "header\sort_tools.h"
# include "header\sgsim.h"

int main() {
    /*
        Parameters:
        (int) mlen = x grid of model
        (int) random seed = random seed
        (int) nR = number of realizations
        (int) hs = lag tolal distance (for variogram calculation)
        (int) bw = steps of lag distance
        (int) vario_flag => if 0 then don't calulate variogram, else 1 then calculate it.
    */

    int mlen = 150;
    int randomseed = 891212;
    int nR = 500;
    int hs = 35;
    int bw = 1;
    int vario_flag = 1;
    srand(randomseed);
    sgsim(mlen, nR, hs, bw, 17.32, 1, vario_flag);

    system("PAUSE");
    return 0;
}
