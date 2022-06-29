# include "../header/sort_tools.hpp"

using namespace std;

template <class T>
void sorting<T>::swap2d(vector<vector<T>>& arr1, int i, int j){

    T data1, data2, data3;

    data1 = arr1[0][j];
    data2 = arr1[1][j];
    data3 = arr1[2][j];

    arr1[0][j] = arr1[0][i];
    arr1[1][j] = arr1[1][i];
    arr1[2][j] = arr1[2][i];

    arr1[0][i] = data1;
    arr1[1][i] = data2;
    arr1[2][i] = data3;
}
template void sorting<double>::swap2d(vector<vector<double>>& arr1, int i, int j);
template void sorting<float>::swap2d(vector<vector<float>>& arr1, int i, int j);
template void sorting<int>::swap2d(vector<vector<int>>& arr1, int i, int j);

template <class T>
int sorting<T>::partition2d(vector<vector<T>>& matrix, int front, int end){

    int mid=front+(end-front)/2;
    
    if(matrix[2][front]>matrix[2][end]){
        
        sorting<T>::swap2d(matrix, front, end);
    }
    if(matrix[2][mid]>matrix[2][end]){
        sorting<T>::swap2d(matrix, mid, end);
    }
    if(matrix[2][mid]>matrix[2][front]){
        sorting<T>::swap2d(matrix, mid, front);
    }
    
    double pivotkey = matrix[2][front];

    while(front < end){

        while(front < end && matrix[2][end] >= pivotkey){
            end--;
        }
        sorting<T>::swap2d(matrix, front, end);

        while(front < end && matrix[2][front] <= pivotkey){
            front++;
        }
        sorting<T>::swap2d(matrix, front, end);
    }
    return front;
}
template int sorting<double>::partition2d(vector<vector<double>>& matrix, int front, int end);
template int sorting<float>::partition2d(vector<vector<float>>& matrix, int front, int end);
template int sorting<int>::partition2d(vector<vector<int>>& matrix, int front, int end);


template <class T>
void sorting<T>::quicksorted_2d(vector<vector<T>>& matrix, int front, int end){

    int pivot;
    
    if(front < end){

        pivot = sorting<T>::partition2d(matrix, front, end);
        sorting<T>::quicksorted_2d(matrix, front, pivot-1);
        sorting<T>::quicksorted_2d(matrix, pivot+1, end);   
    }

}
template void sorting<double>::quicksorted_2d(vector<vector<double>>& matrix, int front, int end);
template void sorting<float>::quicksorted_2d(vector<vector<float>>& matrix, int front, int end);
template void sorting<int>::quicksorted_2d(vector<vector<int>>& matrix, int front, int end);

