/**
 * @file sgsim.h
 * @brief Header file for Sequential Gaussian Simulation (SGSIM) library.
 *
 * This header defines the data structure sgsim_t and several functions for
 * conducting sequential Gaussian simulation.
 *
 * Copyright (c) 2022-2023 Zncl2222
 * License: MIT
 */

#ifndef UC_SGSIM_C_CORE_INCLUDE_SGSIM_H_
#define UC_SGSIM_C_CORE_INCLUDE_SGSIM_H_

# include "cov_model.h"
# include "../c_array_tools/src/c_array.h"

/**
 * @struct sgsim_t
 * @brief Structure to hold parameters and data for SGSIM simulation.
 *
 * This structure contains various parameters and data needed for performing
 * Sequential Gaussian Simulation (SGSIM).
 */
typedef struct {
    int x_len;  // Length of the realization in the x-direction.
    int realization_numbers;  // Number of realizetions to generate.
    int randomseed;
    int kriging_method;
    int if_alloc_memory;  // Flag to indicate memory allocation status
    int iteration_limit;  // The tolerance of maximum times of iteration error
    double* array;  // Array to store the simulated values.
    double z_min;  // Minimum simulated value.
    double z_max;  // Maximum simulated value.
} sgsim_t;

/**
 * @brief Set default values for an sgsim_t structure.
 *
 * @param sgsim Pointer to an sgsim_t structure to set defaults for.
 * @param cov_model Pointer to the covariance model to be used.
 */
void set_sgsim_defaults(sgsim_t* sgsim, cov_model_t* cov_model);

/**
 * @brief Run Sequential Gaussian Simulation (SGSIM).
 *
 * This function performs Sequential Gaussian Simulation (SGSIM) based on the
 * provided sgsim_t parameters and covariance model.
 *
 * @param sgsim Pointer to the sgsim_t structure containing simulation parameters.
 * @param cov_model Pointer to the covariance model to be used.
 * @param vario_flag Flag indicating whether to calculate variogram or not.
 */
void sgsim_run(sgsim_t* sgsim, cov_model_t* cov_model, int vario_flag);

/**
 * @brief Free memory allocated for an sgsim_t structure.
 *
 * This function releases memory allocated for the sgsim_t structure.
 *
 * @param sgsim Pointer to the sgsim_t structure to free.
 */
void sgsim_t_free(sgsim_t* sgsim);

/**
 * @brief Internal function to free memory used by the SGSIM library.
 *
 * This function is for internal use and frees any additional memory used
 * by the SGSIM library. Users typically don't need to call this directly.
 */
static void sgsim_memory_free();

#endif   // UC_SGSIM_C_CORE_INCLUDE_SGSIM_H_
