// Copyright 2022 Zncl2222
#ifndef UC_SGSIM_C_CORE_INCLUDE_KRIGE_H_
#define UC_SGSIM_C_CORE_INCLUDE_KRIGE_H_

# include "../lib/c_array.h"

void Krige_paramsetting(double a, double C0);

void SimpleKrige(double* array, double* sampled, double* u_array,
                int array_size, double unsampled_point, int idx,
                int neighbor, mt19937_state* rng_state);

int find_neighbor(double* array, double* sampled, double* u_array,
                int array_size, double unsampled_point, int idx,
                int neighbor, mt19937_state* rng_state);

void krige_memory_free();

#endif  // UC_SGSIM_C_CORE_INCLUDE_KRIGE_H_
