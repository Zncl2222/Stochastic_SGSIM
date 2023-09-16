// Copyright 2022 Zncl2222
#ifndef UC_SGSIM_C_CORE_INCLUDE_SORT_TOOLS_H_
#define UC_SGSIM_C_CORE_INCLUDE_SORT_TOOLS_H_

void swaprows(double** x, int i, int j);

void quicksort2d(double** array, int front, int end);

int partition2d(double** array, int front, int end);

void quickselect2d(double** array, int front, int end, int k);

#endif  // UC_SGSIM_C_CORE_INCLUDE_SORT_TOOLS_H_
