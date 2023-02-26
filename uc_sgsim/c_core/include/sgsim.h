// Copyright 2022 Zncl2222
#ifndef UC_SGSIM_C_CORE_INCLUDE_SGSIM_H_
#define UC_SGSIM_C_CORE_INCLUDE_SGSIM_H_

void sgsim(int X, int nR, int hs, int bw,
        double range, double sill, int randomseed, int vario_flag);

void sgsim_dll(double* RandomFieldX, int X, int nR, double a, double C0, int randomseed);

void sgsim_memory_free();

#endif   // UC_SGSIM_C_CORE_INCLUDE_SGSIM_H_
