// Copyright 2022 Zncl2222
#ifndef UC_SGSIM_C_CORE_INCLUDE_COV_MODEL_H_
#define UC_SGSIM_C_CORE_INCLUDE_COV_MODEL_H_

typedef struct {
    int bw;
    int hs;
    double range;
    double sill;
} cov_model_t;

void cov_model_init(
    cov_model_t* _cov_model, int bw,
    int hs, double range, double sill);

void cov_model(const double *x, double* cov, int n_dim, double a, double C0);

void cov_model2d(const double **x, double* cov, int n_dim, double a, double C0);

#endif  // UC_SGSIM_C_CORE_INCLUDE_COV_MODEL_H_
