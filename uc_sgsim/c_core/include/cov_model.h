// Copyright 2022 Zncl2222
#ifndef UC_SGSIM_C_CORE_INCLUDE_COV_MODEL_H_
#define UC_SGSIM_C_CORE_INCLUDE_COV_MODEL_H_

typedef struct {
    int bw_l;
    int bw_s;
    int bw;
    double k_range;
    double sill;
    double nugget;
} cov_model_t;

void cov_model_init(
    cov_model_t* _cov_model, int bw_l,
    int bw_s, double k_range, double sill, double nugget);

void cov_model(const double* x, double* cov, int n_dim, cov_model_t* _cov_model);

void cov_model2d(const double** x, double* cov, int n_dim, cov_model_t* _cov_model);

#endif  // UC_SGSIM_C_CORE_INCLUDE_COV_MODEL_H_
