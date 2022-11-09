// Copyright 2022 Zncl2222
#ifndef UC_SGSIM_C_CORE_HEADER_SORT_TOOLS_H_
#define UC_SGSIM_C_CORE_HEADER_SORT_TOOLS_H_
#include <utility>

double** sort2d(double** x, int n_dim);

void swap(double** x, int i, int j);

double** BubbleSort2d(double** x, int n_dim);

void quicksort2d(double** array, int front, int end);

int Partition2d(double** array, int front, int end);

void quicksort1d(double* array, int front, int end);

int Partition1d(double* array, int front, int end);

void swap1d(double* X, double* Y);

#endif  // UC_SGSIM_C_CORE_HEADER_SORT_TOOLS_H_
