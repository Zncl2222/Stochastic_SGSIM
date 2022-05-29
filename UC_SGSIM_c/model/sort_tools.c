# include "..\header\sort_tools.h"
# include <stdio.h>


void swap(double** x, int i, int j){

    double tempp;
    double tempp2;
    double tempp3;

    tempp = x[j][0];
    tempp2= x[j][1];
    tempp3= x[j][2];

    x[j][0] = x[i][0];
    x[j][1] = x[i][1];
    x[j][2] = x[i][2];

    x[i][0] = tempp;
    x[i][1] = tempp2;
    x[i][2] = tempp3;

}

void swap1d(double* X, double* Y){

    double temp = *X;
    *X = *Y;
    *Y = temp;
}

double** sort2d(double** x, int n_dim){   

    for(int i = 0; i < n_dim; i++) {

       for(int j = i + 1; j < n_dim; j++) {

           if( x[j][2] < x[i][2] ) {

               swap(x, i , j);
           }
       }
   }
   return x;
}

double** BubbleSort2d(double** x, int n_dim){

    for(int i = 0; i < n_dim; i ++){

        for(int j = n_dim - 1; j >= i; j--){

            if(x[j][2] > x[j+1][2]){
                swap(x, j, j+1);
            }
        }

    }
    return x;
}

int Partition2d(double** array, int front, int end){

    int mid=front+(end-front)/2;
        
    if(array[front][2]>array[end][2]){
        swap(array, front, end);
    }
    if(array[mid][2]>array[end][2]){
        swap(array, mid, end);
    }
    if(array[mid][2]>array[front][2]){
        swap(array, mid, front);
    }
    
    double pivotkey = array[front][2];

    while(front < end){

        while(front < end && array[end][2] >= pivotkey){
            end--;
        }
        swap(array, front, end);

        while(front < end && array[front][2] <= pivotkey){
            front++;
        }
        swap(array, front, end);
    }
    return front;
}

void quicksort2d(double** array, int front, int end){

    int pivot;

    if(front < end){

        pivot = Partition2d(array, front, end);
        quicksort2d(array, front, pivot-1);
        quicksort2d(array, pivot+1, end);   
    }
}

int Partition1d(double* array, int front, int end){

    double pivotkey = array[front];

    while(front < end){

        while(front < end && array[end] >= pivotkey){
            end--;
        }
        swap1d(&array[front], &array[end]);

        while(front < end && array[front] <= pivotkey){
            front++;
        }
        swap1d(&array[front], &array[end]);
    }
    return front;
}


void quicksort1d(double* array, int front, int end){

    int pivot;

    if(front < end){

        pivot = Partition1d(array, front, end);
        quicksort1d(array, front, pivot-1);
        quicksort1d(array, pivot+1, end);   
    }
}

