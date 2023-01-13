// Copyright 2022 Zncl2222
#ifndef UC_SGSIM_C_CORE_HEADER_SGSIM_H_
#define UC_SGSIM_C_CORE_HEADER_SGSIM_H_

void sgsim(int X, int nR, int hs, int bw,
        double range, double sill, int vario_flag);

void sgsim_dll(double* RandomFieldX, int X, int nR, double a, double C0);

void sgsim_memory_free();

#endif   // UC_SGSIM_C_CORE_HEADER_SGSIM_H_
