is_excute=false
file_name="uc_sgsim"

if $is_excute; then
    # compile executable
    echo "Building executable file"
    file=$file_name".out"
    gcc -o uc_sgsim.out \
    ../../UC_SGSIM_c_example.c \
    src/sgsim.c  src/kriging.c  src/cov_model.c \
    src/variogram.c  src/sort_tools.c  src/random_tools.c  src/matrix_tools.c -lm
else
    # compile shared lib
    echo "Building shared library"
    file=$file_name".so"
    gcc -shared -fPIC -o $file  \
    src/sgsim.c  src/kriging.c  src/cov_model.c \
    src/variogram.c  src/sort_tools.c  src/random_tools.c  src/matrix_tools.c
fi
