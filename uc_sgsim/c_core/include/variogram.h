// Copyright 2022 Zncl2222
#ifndef UC_SGSIM_C_CORE_INCLUDE_VARIOGRAM_H_
#define UC_SGSIM_C_CORE_INCLUDE_VARIOGRAM_H_

void variogram(const double* array, double* v, int mlen, int bw, int bw_s);

double variance(const double* array, int mlen);

#endif  // UC_SGSIM_C_CORE_INCLUDE_VARIOGRAM_H_
