// Copyright 2022 Zncl2222
#ifndef UC_SGSIM_C_CORE_HEADER_COV_MODEL_H_
#define UC_SGSIM_C_CORE_HEADER_COV_MODEL_H_

void Cov_model(double *x, double* cov, int n_dim, double a, double C0);

void Cov_model2d(double **x, double* cov, int n_dim, double a, double C0);

#endif  // UC_SGSIM_C_CORE_HEADER_COV_MODEL_H_
