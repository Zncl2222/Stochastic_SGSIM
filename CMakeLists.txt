cmake_minimum_required(VERSION 2.4.0)
enable_testing()

project(uc_sgsim_c_example)
set(CMAKE_VERBOSE_MAKEFILE ON)
# set true to compile the executable file (win -> .exe, linux -> .out)
set(IS_EXECUTE true)

# set output directory
if(IS_EXECUTE)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
else()
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/uc_sgsim/c_core")
endif()

include_directories(${CMAKE_SOURCE_DIR})

# run CMakeLists in uc_sgsim folder
add_subdirectory(uc_sgsim)

# compile and install executable file
if(IS_EXECUTE)
    add_executable(uc_sgsim_c_example UC_SGSIM_c_example.c)
    target_link_libraries(uc_sgsim_c_example uc_sgsim)
    if(WIN32)
        set_target_properties(uc_sgsim_c_example PROPERTIES SUFFIX ".exe")
    else()
        set_target_properties(uc_sgsim_c_example PROPERTIES SUFFIX ".out")
    endif()

    # make sure example can run
    add_test(
        uc_sgsim_c_example
        ${CMAKE_COMMAND}
        --build ${CMAKE_BINARY_DIR}
        --target uc_sgsim_c_example
    )
    # memory check (valgrind)
    find_program(VALGRIND "valgrind")

    if(VALGRIND)
        add_custom_target(
            memcheck
            COMMAND ${VALGRIND}
            --tool=memcheck
            --leak-check=full
            --track-origins=yes
            --show-reachable=yes
            --error-exitcode=1
            $<TARGET_FILE:uc_sgsim_c_example>
        )
    endif()
endif()
