gcc -c model/sgsim.c
gcc -c model/krige.c
gcc -c model/cov_model.c
gcc -c model/variogram.c
gcc -c model/sort_tools.c
gcc -c model/random_tools.c
gcc -c model/matrix_tools.c
gcc -c UC_SGSIM_c_example.c

gcc sgsim.o krige.o cov_model.o variogram.o sort_tools.o random_tools.o matrix_tools.o UC_SGSIM_c_example.o -o uc_sgsim.exe
