// Copyright 2022 Zncl2222

# include <stdio.h>

# include "../include/sort_tools.h"


void swaprows(double** x, int i, int j) {
    double* temp = x[j];
    x[j] = x[i];
    x[i] = temp;
}

int partition2d(double** array, int front, int end) {
    double pivotkey = array[front][2];
    int i = front - 1;
    int j = end + 1;

    while (1) {
        do {
            i++;
        } while (array[i][2] < pivotkey);

        do {
            j--;
        } while (array[j][2] > pivotkey);

        if (i >= j) {
            return j;
        }

        swaprows(array, i, j);
    }
}

void quicksort2d(double** array, int front, int end) {
    if (front < end) {
        int pivot = partition2d(array, front, end);
        quicksort2d(array, front, pivot);
        quicksort2d(array, pivot + 1, end);
    }
}

void quickselect2d(double** array, int front, int end, int k) {
    if (front < end) {
        int pivot = partition2d(array, front, end);

        if (pivot == k - 1) {
            return;
        } else if (pivot > k - 1) {
            quickselect2d(array, front, pivot, k);
        } else {
            quickselect2d(array, pivot + 1, end, k);
        }
    }
}
