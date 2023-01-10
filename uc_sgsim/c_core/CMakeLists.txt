cmake_minimum_required(VERSION 2.4.0)

project(c_core)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/uc_sgsim/c_core/")

if(NOT DEFINED IS_EXECUTE OR NOT IS_EXECUTE)
    set(LIB_TYPE SHARED)
else()
    set(LIB_TYPE STATIC)
endif()

# Include header files
include_directories(header/)

# Find source files
add_library(
    uc_sgsim ${LIB_TYPE}
    model/cov_model.c
    model/krige.c
    model/matrix_tools.c
    model/random_tools.c
    model/sgsim.c
    model/sort_tools.c
    model/variogram.c
)

if(NOT DEFINED IS_EXECUTE OR NOT IS_EXECUTE)
    install(TARGETS uc_sgsim DESTINATION "${CMAKE_SOURCE_DIR}/uc_sgsim/c_core/")
    set_target_properties(uc_sgsim PROPERTIES PREFIX "")
else()
    target_link_libraries(uc_sgsim PRIVATE m)
endif()
