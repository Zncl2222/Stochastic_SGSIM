/**
 * @file kriging.h
 * @brief Header file for Kriging functions in the SGSIM library.
 *
 * This header defines the data structures and functions related to Kriging
 * used in the Sequential Gaussian Simulation (SGSIM) library.
 *
 * Copyright (c) 2022-2023 Zncl2222
 * License: MIT
 */

#ifndef UC_SGSIM_C_CORE_INCLUDE_KRIGING_H_
#define UC_SGSIM_C_CORE_INCLUDE_KRIGING_H_

# include "../lib/c_array.h"
# include "../include/cov_model.h"

/**
 * @struct sampling_state
 * @brief Structure to represent the state of sampling during Kriging.
 *
 * This structure maintains the state of sampling during the Kriging process,
 * including information about neighbors, current length, and sampled points.
 */
typedef struct {
    int neighbor;             // Number of neighbors considered.
    int currlen;              // Current length of sampled points.
    int idx;                  // Index used for sampling state.
    double unsampled_point;   // Unsampling point during Kriging.
    c_array_double sampled;   // Array to store sampled points.
    c_array_double u_array;   // Array to store unsampled points.
} sampling_state;

/**
 * @brief Initialize a sampling state structure.
 *
 * This function initializes a sampling_state structure with default values.
 *
 * @param sampling Pointer to the sampling_state structure to initialize.
 * @param x_grid_len Length of the grid used in sampling.
 */
void sampling_state_init(sampling_state* sampling, int x_grid_len);

/**
 * @brief Update the sampling state with a new unsampled point.
 *
 * This function updates the sampling state with a new unsampled point and
 * its index in the grid.
 *
 * @param sampling Pointer to the sampling_state structure to update.
 * @param unsampled_point The unsampled point to add.
 * @param idx Index of the unsampled point in the grid.
 */
void sampling_state_update(sampling_state* sampling, double unsampled_point, int idx);

/**
 * @brief Set Kriging parameters for the simulation.
 *
 * This function sets the Kriging parameters based on the length of the grid
 * and the covariance model to be used.
 *
 * @param x_len Length of the grid.
 * @param cov_model Pointer to the covariance model.
 */
void kriging_param_setting(int x_len, const cov_model_t* cov_model);

/**
 * @brief Perform simple Kriging to estimate values at unsampled points.
 *
 * This function performs simple Kriging to estimate values at unsampled points
 * based on the sampling state, random number generator state, and Kriging method.
 *
 * @param array Array to store estimated values.
 * @param sampling Pointer to the sampling state.
 * @param rng_state Pointer to the random number generator state.
 * @param kriging_method Kriging method to use (e.g., ordinary kriging).
 */
void simple_kriging(double* array, sampling_state* sampling,
                    mt19937_state* rng_state, int kriging_method,
                    int use_cov_cache);

/**
 * @brief Find neighbor points for Kriging.
 *
 * This function finds neighbor points for Kriging estimation based on the
 * sampling state.
 *
 * @param array Array of values.
 * @param sampling Pointer to the sampling state.
 * @param rng_state Pointer to the random number generator state.
 * @return true or false by 1 or 0 (true means there is a neighbor for kriging to interpolate)
 */
int find_neighbor(double* array, sampling_state* sampling, mt19937_state* rng_state);

/**
 * @brief Augment a matrix for Ordinary Kriging.
 *
 * This function augments a matrix for Ordinary Kriging estimation based on the number
 * of neighbor points.
 *
 * @param mat Matrix to be augmented.
 * @param neighbor Number of neighbor points.
 */
void matrix_augmented(double** mat, int neighbor);

/**
 * @brief Free memory used for Kriging.
 *
 * This function frees memory used by the Kriging functions in the SGSIM library.
 */
void kriging_memory_free();

#endif  // UC_SGSIM_C_CORE_INCLUDE_KRIGING_H_
