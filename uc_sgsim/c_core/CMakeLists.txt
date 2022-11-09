cmake_minimum_required(VERSION 2.4.0)

project(UC_SGSIM_c)

# Find source files
add_library(
    UC_SGSIM_c SHARED
    UC_SGSIM_c_example.c
    model/cov_model.c
    model/krige.c
    model/matrix_tools.c
    model/random_tools.c
    model/sgsim.c
    model/sort_tools.c
    model/variogram.c
)

# Include header files
include_directories(header/)

# Install library
install(TARGETS UC_SGSIM_c DESTINATION lib/UC_SGSIM_c)

# Install library headers
add_library(
    UC_SGSIM_c SHARED
    header/cov_model.h
    header/krige.h
    header/matrix_tools.h
    header/random_tools.h
    header/sgsim.h
    header/sort_tools.h
    header/variogram.h
)
install(FILES ${HEADERS} DESTINATION include/UC_SGSIM_c)
