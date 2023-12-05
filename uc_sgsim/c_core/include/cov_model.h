/**
 * @file cov_model.h
 * @brief Declarations for Covariance Model Functions
 *
 * This file contains declarations for functions related to covariance models in geostatistics.
 * These functions are used to set default values and calculate covariance values based on
 * a specified model.
 *
 * Copyright (c) 2022-2023 Zncl2222
 * License: MIT
 */

#ifndef UC_SGSIM_C_CORE_INCLUDE_COV_MODEL_H_
#define UC_SGSIM_C_CORE_INCLUDE_COV_MODEL_H_

/**
 * @struct cov_model_t
 * @brief Structure to hold covariance model parameters.
 *
 * This structure defines the parameters of a covariance model used in geostatistical calculations.
 */
typedef struct {
    int bw_l;              // bandwidth length
    int bw_s;              // bandwidth steps
    int bw;                // bandwidth
    int max_neighbor;      // maximum number of neighbors
    int use_cov_cache;     // flag to indicate whether to use covariance cache
    double k_range;        // kriging range
    double sill;           // sill value
    double nugget;         // nugget value
} cov_model_t;

/**
 * @brief Set default values for a covariance model.
 *
 * This function sets default values for a covariance model if they are not already specified.
 *
 * @param _cov_model A pointer to a cov_model_t structure to be initialized.
 */
void set_cov_model_default(cov_model_t* _cov_model);

/**
 * @brief Calculate covariance for a one-dimensional dataset.
 *
 * This function calculates the covariance values for a one-dimensional dataset based on the
 * specified covariance model parameters.
 *
 * @param x An array of double values representing the spatial coordinates.
 * @param cov An array to store the calculated covariance values.
 *            The size should be at least n_dim (Results will be save at here).
 * @param n_dim The number of dimensions (length of the x array).
 * @param _cov_model A pointer to a cov_model_t structure that defines the covariance model parameters.
 */
void cov_model(const double* x, double* cov, int n_dim, cov_model_t* _cov_model);

/**
 * @brief Calculate covariance for a two-dimensional dataset.
 *
 * This function calculates the covariance values for a two-dimensional dataset based on the
 * specified covariance model parameters.
 *
 * @param x A 2D array of double values representing the spatial coordinates. The dimensions
 * should be n_dim x n_dim.
 * @param cov An array to store the calculated covariance values.
 *            The size should be at least n_dim * n_dim (Results will be save at here).
 * @param n_dim The number of dimensions (length of each side of the 2D array).
 * @param _cov_model A pointer to a cov_model_t structure that defines the covariance model parameters.
 */
void cov_model2d(const double** x, double* cov, int n_dim, cov_model_t* _cov_model);

#endif  // UC_SGSIM_C_CORE_INCLUDE_COV_MODEL_H_
