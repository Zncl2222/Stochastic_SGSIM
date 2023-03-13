// Copyright 2022 Zncl2222

# include <stdio.h>

# include "../include/sort_tools.h"


void swap(double** x, int i, int j) {
    double tempp;
    double tempp2;
    double tempp3;

    tempp = x[j][0];
    tempp2 = x[j][1];
    tempp3 = x[j][2];

    x[j][0] = x[i][0];
    x[j][1] = x[i][1];
    x[j][2] = x[i][2];

    x[i][0] = tempp;
    x[i][1] = tempp2;
    x[i][2] = tempp3;
}

int Partition2d(double** array, int front, int end) {
    int mid = front + (end-front) / 2;

    if (array[front][2] > array[end][2]) {
        swap(array, front, end);
    }
    if (array[mid][2] > array[end][2]) {
        swap(array, mid, end);
    }
    if (array[mid][2] > array[front][2]) {
        swap(array, mid, front);
    }

    double pivotkey = array[front][2];

    while (front < end) {
        while (front < end && array[end][2] >= pivotkey) {
            end--;
        }
        swap(array, front, end);

        while (front < end && array[front][2] <= pivotkey) {
            front++;
        }
        swap(array, front, end);
    }
    return front;
}

void quicksort2d(double** array, int front, int end) {
    int pivot;

    if (front < end) {
        pivot = Partition2d(array, front, end);
        quicksort2d(array, front, pivot-1);
        quicksort2d(array, pivot+1, end);
    }
}
