// Copyright 2022 Zncl2222
#ifndef UC_SGSIM_C_CORE_INCLUDE_SGSIM_H_
#define UC_SGSIM_C_CORE_INCLUDE_SGSIM_H_

# include "cov_model.h"
# include "../lib/c_array.h"

typedef struct {
    int x_len;
    int realization_numbers;
    int randomseed;
    int if_alloc_memory;
    double* array;
} sgsim_t;

void sgsim_init(sgsim_t* _sgsim, int x_grid,
                int realization_numbers, int randomseed, int if_alloc_memory);

void sgsim_run(sgsim_t* _sgsim, const cov_model_t* cov_model, int vario_flag);

void sgsim_t_free(sgsim_t* _sgsim);

static void sgsim_memory_free();

#endif   // UC_SGSIM_C_CORE_INCLUDE_SGSIM_H_
