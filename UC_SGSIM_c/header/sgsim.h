#ifndef SGSIM_H
#define SGSIM_H

void sgsim(int X, int nR, int hs, int bw, double range, double sill, int randomseed, int vario_flag);

void sgsim_dll(double* RandomFieldX, int X, int nR, double a, double C0, int randomseed);

void free_memory(int mlen);

#endif