cmake_minimum_required(VERSION 2.4.0)

project(uc_sgsim_c_example)

set(IS_EXECUTE FALSE)

if(IS_EXECUTE)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
else()
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/uc_sgsim/c_core")
endif()

include_directories(${CMAKE_SOURCE_DIR})

add_subdirectory(uc_sgsim)

if(IS_EXECUTE)
add_executable(uc_sgsim_c_example UC_SGSIM_c_example.c)

target_link_libraries(uc_sgsim_c_example uc_sgsim)

set_target_properties(uc_sgsim_c_example PROPERTIES SUFFIX ".out")
endif()
