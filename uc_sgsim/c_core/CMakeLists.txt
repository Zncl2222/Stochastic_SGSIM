cmake_minimum_required(VERSION 2.4.0)

project(uc_sgsim)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Find source files
add_library(
    uc_sgsim SHARED
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
install(TARGETS uc_sgsim)

set_target_properties(uc_sgsim PROPERTIES PREFIX "")
