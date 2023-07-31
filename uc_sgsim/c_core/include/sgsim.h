// Copyright 2022 Zncl2222
#ifndef UC_SGSIM_C_CORE_INCLUDE_SGSIM_H_
#define UC_SGSIM_C_CORE_INCLUDE_SGSIM_H_

# include "cov_model.h"
# include "../lib/c_array.h"

typedef struct {
    int x_len;
    int realization_numbers;
    int randomseed;
    int kriging_method;
    int if_alloc_memory;
    double* array;
    double z_min;
    double z_max;
} sgsim_t;


void set_sgsim_defaults(sgsim_t* _sgsim, cov_model_t* cov_model);

void sgsim_run(sgsim_t* _sgsim, cov_model_t* cov_model, int vario_flag);

void sgsim_t_free(sgsim_t* _sgsim);

static void sgsim_memory_free();

#endif   // UC_SGSIM_C_CORE_INCLUDE_SGSIM_H_
