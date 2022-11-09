// Copyright 2022 Zncl2222

# include <stdio.h>
# include <stdlib.h>

# include ".\uc_sgsim\c_core\header\matrix_tools.h"
# include ".\uc_sgsim\c_core\header\variogram.h"
# include ".\uc_sgsim\c_core\header\krige.h"
# include ".\uc_sgsim\c_core\header\cov_model.h"
# include ".\uc_sgsim\c_core\header\sort_tools.h"
# include ".\uc_sgsim\c_core\header\sgsim.h"

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

    sgsim(mlen, nR, hs, bw, 17.32, 1, randomseed, vario_flag);

    system("PAUSE");
    return 0;
}
