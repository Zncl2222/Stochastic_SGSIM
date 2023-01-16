// Copyright 2022 Zncl2222
#ifndef UC_SGSIM_C_CORE_INCLUDE_VARIOGRAM_H_
#define UC_SGSIM_C_CORE_INCLUDE_VARIOGRAM_H_

void variogram(double* array, double* v, int mlen, int hs, int steps);

double variance(double* array, int mlen);

#endif  // UC_SGSIM_C_CORE_INCLUDE_VARIOGRAM_H_
