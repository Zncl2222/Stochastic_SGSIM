/**
 * @file sort_tools.h
 * @brief Sorting Tools Library
 *
 * This header file provides declarations for functions related to sorting and manipulation of
 * two-dimensional arrays. The functions in this library assist in sorting rows of a two-dimensional
 * array and selecting elements using the quicksort and quickselect algorithms.
 *
 * Copyright (c) 2022 Zncl2222
 */

#ifndef UC_SGSIM_C_CORE_INCLUDE_SORT_TOOLS_H_
#define UC_SGSIM_C_CORE_INCLUDE_SORT_TOOLS_H_

/**
 * @brief Swap two rows in a two-dimensional array.
 *
 * This function swaps two rows in a two-dimensional array. It is commonly used in sorting algorithms.
 *
 * @param x A pointer to a two-dimensional array of double values.
 * @param i The index of the first row to swap.
 * @param j The index of the second row to swap.
 */
void swaprows(double** x, int i, int j);

/**
 * @brief Sort a two-dimensional array using the quicksort algorithm.
 *
 * This function sorts a two-dimensional array based on the values in a specified column using the
 * quicksort algorithm.
 *
 * @param array A pointer to a two-dimensional array of double values.
 * @param front The index of the first row to start sorting.
 * @param end The index of the last row to end sorting.
 */
void quicksort2d(double** array, int front, int end);

/**
 * @brief Partition a two-dimensional array for the quicksort algorithm.
 *
 * This function partitions a two-dimensional array based on the values in a specified column as
 * part of the quicksort algorithm.
 *
 * @param array A pointer to a two-dimensional array of double values.
 * @param front The index of the first row to start partitioning.
 * @param end The index of the last row to end partitioning.
 *
 * @return The index of the pivot element.
 */
int partition2d(double** array, int front, int end);

/**
 * @brief Select an element from a two-dimensional array using the quickselect algorithm.
 *
 * This function selects an element from a two-dimensional array based on the values in a specified column using the
 * quickselect algorithm.
 *
 * @param array A pointer to a two-dimensional array of double values.
 * @param front The index of the first row to start sorting.
 * @param end The index of the last row to end sorting.
 * @param k The number of the element to select.
 */
void quickselect2d(double** array, int front, int end, int k);

#endif  // UC_SGSIM_C_CORE_INCLUDE_SORT_TOOLS_H_
