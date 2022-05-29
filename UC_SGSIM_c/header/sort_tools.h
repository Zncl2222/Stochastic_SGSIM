#ifndef SORT_TOOLS_H
#define SORT_TOOLS_H

double** sort2d(double** x, int n_dim);

void swap(double** x, int i, int j);

double** BubbleSort2d(double** x, int n_dim);

void quicksort2d(double** array, int front, int end);

int Partition2d(double** array, int front, int end);

void quicksort1d(double* array, int front, int end);

int Partition1d(double* array, int front, int end);

void swap1d(double* X, double* Y);

#endif