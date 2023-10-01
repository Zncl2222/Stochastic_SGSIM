/**
 * @file variogram.h
 * @brief Variogram Calculation Library
 *
 * This header file provides declarations for functions related to variogram calculation.
 * The functions in this library assist in calculating variograms and variances from input data.
 *
 * Copyright (c) 2022-2023 Zncl2222
 * License: MIT
 */
#ifndef UC_SGSIM_C_CORE_INCLUDE_VARIOGRAM_H_
#define UC_SGSIM_C_CORE_INCLUDE_VARIOGRAM_H_

/**
 * @brief Calculate the experimental variogram.
 *
 * This function calculates the experimental variogram of a given data array using specified bandwidths.
 *
 * @param array A pointer to the input data array.
 * @param v A pointer to the output variogram array.
 * @param mlen The length of the input data array.
 * @param bw The bandwidth used for the variogram calculation.
 * @param bw_s The bandwidth step used for the variogram calculation.
 */
void variogram(const double* array, double* v, int mlen, int bw, int bw_s);

/**
 * @brief Calculate the variance of a data array.
 *
 * This function calculates the variance of a given data array.
 *
 * @param array A pointer to the input data array.
 * @param mlen The length of the input data array.
 *
 * @return The calculated variance of the input data.
 */
double variance(const double* array, int mlen);

#endif  // UC_SGSIM_C_CORE_INCLUDE_VARIOGRAM_H_
