// Copyright 2022 Zncl2222
#ifndef UC_SGSIM_C_CORE_INCLUDE_SGSIM_H_
#define UC_SGSIM_C_CORE_INCLUDE_SGSIM_H_

# include "cov_model.h"
# include "../lib/c_array.h"

struct sgsim_t {
    int x_grid;
    int realization_numbers;
    int randomseed;
    int if_alloc_memory;
    double* array;
};

void sgsim_init(struct sgsim_t* _sgsim, int x_grid,
                int realization_numbers, int randomseed, int if_alloc_memory);

void sgsim_run(struct sgsim_t* _sgsim, const struct cov_model_t* cov_model, int vario_flag);

void sgsim_t_free(struct sgsim_t* _sgsim);

static void sgsim_memory_free();

#endif   // UC_SGSIM_C_CORE_INCLUDE_SGSIM_H_
