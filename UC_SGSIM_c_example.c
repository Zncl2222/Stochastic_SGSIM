# include <stdio.h>
# include <stdlib.h>

# include ".\UC_SGSIM_c\header\matrix_tools.h"
# include ".\UC_SGSIM_c\header\variogram.h"
# include ".\UC_SGSIM_c\header\krige.h"
# include ".\UC_SGSIM_c\header\cov_model.h"
# include ".\UC_SGSIM_c\header\sort_tools.h"
# include ".\UC_SGSIM_c\header\sgsim.h"

/* Developer's set : mingw32 with gcc version 9.2.0*/
/* Compile each *.c file to object(*.o) first (don't compile the header file (*.h))*/
/* Then compile all .o file become the .exe*, for instance : */
/* "gcc -c sgsim.c","gcc -c krige.c", "gcc -c cov_model.c" etc........*/
/* "gcc Uc_SGSIM3.o sgsim.o krige.o cov_model.o matrix_tools.o random_tools.o sort_tools.o variogram.o -o name.exe"*/

int main(){

    /*  Parameters:
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