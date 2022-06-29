
# include "../header/matrix_tools.hpp"
# include <vector>
# include <cstdio>
# include <string.h>
//# include <sys/io.h>  // for linux
//# include <sys/stat.h> // for linux
# include <io.h>     // for windows

using namespace std;

template<class T>
void matrix_tools<T>::LUdecomposition(vector<vector<T>> a ,vector<T> b, vector<T>& x , int n){

   int i = 0, j = 0, k = 0;
   vector<vector<T>> l(10,vector<T>(10));
   vector<vector<T>> u(10,vector<T>(10));
   vector<T> c(10);

   for (i = 0; i < n; i++) {
       
        for (j = 0; j < n; j++) {

            if (j < i)
                l[j][i] = 0;
            else {

                l[j][i] = a[j][i];
                for (k = 0; k < i; k++) {
                    l[j][i] = l[j][i] - l[j][k] * u[k][i];
                }
            }
        }
      for (j = 0; j < n; j++) {

        if (j < i)
            u[i][j] = 0;
        else if (j == i)
            u[i][j] = 1;
        else{
            u[i][j] = a[i][j] / l[i][i];
            for (k = 0; k < i; k++) {
               u[i][j] = u[i][j] - ((l[i][k] * u[k][j]) / l[i][i]);
            }
         }
      }
   }
    // Solve L(Ux)=b, assume Ux=c 
    c[0] = b[0] / l[0][0];
    
    for (i = 1; i < n; i++){

        c[i] = b[i];
        
        for (j = 0; j < i; j++){

            c[i] = c[i] - l[i][j] * c[j]; 
        }

        c[i] = c[i] / l[i][i];
    }

    x[n] = c[n];

    for (i = n - 1; i >= 0; i--){

        x[i] = c[i];

        for (j = i + 1; j < n; j++){
            x[i] -= u[i][j] * x[j];
        }
    }

}
template void matrix_tools<double>::
LUdecomposition(vector<vector<double>> a ,vector<double> b, vector<double>& x , int n);
template void matrix_tools<float>::
LUdecomposition(vector<vector<float>> a ,vector<float> b, vector<float>& x , int n);
template void matrix_tools<int>::
LUdecomposition(vector<vector<int>> a ,vector<int> b, vector<int>& x , int n);


template<class T>
vector<T> matrix_tools<T>::arange(int x){
    
    vector<T> res(x);
    //this->array_ = vector<T>(x);

    for (int i = 0; i < x; i++){
        res[i] = i;
    }
    
    return res;
}
template vector<double> matrix_tools<double>::arange(int x);
template vector<float> matrix_tools<float>::arange(int x);
template vector<int> matrix_tools<int>::arange(int x);

template<class T>
void matrix_tools<T>::pdist(vector<T> x, vector<vector<T>>& c, int n_dim){

    if(n_dim < 1)
        return;

    for (int i = 0; i < n_dim; i++){

        for (int j = 0; j < n_dim; j++){

            c[i][j] = fabs(x[j] - x[i]);
        }
    }

}
template void matrix_tools<double>::pdist(vector<double> x, vector<vector<double>>& c, int n_dim);
template void matrix_tools<float>::pdist(vector<float> x, vector<vector<float>>& c, int n_dim);
template void matrix_tools<int>::pdist(vector<int> x, vector<vector<int>>& c, int n_dim);

template<class T>
void matrix_tools<T>::matrixform(vector<T> x,vector<vector<T>>& matrix, int n_dim){

 
    for (int i = 0; i < n_dim; i++){

        for (int j = 0; j < n_dim; j++){

            matrix[i][j] = x[n_dim*i+j];
        }
    }
}
template void matrix_tools<double>::matrixform(vector<double> x,vector<vector<double>>& matrix, int n_dim);
template void matrix_tools<float>::matrixform(vector<float> x,vector<vector<float>>& matrix, int n_dim);
template void matrix_tools<int>::matrixform(vector<int> x,vector<vector<int>>& matrix, int n_dim);

template<class T>
vector<vector<T>> matrix_tools<T>::matrixReshape(vector<vector<T>> mat, int r, int c){
    
    int matSize = mat.size();
    int matCol = mat[0].size();

    if(matSize*matCol != r*c || (matSize==r && matCol==c)){
        printf("(matrixReshape) ERROR: Matrix can't reshape to the target shape.\n");
        return mat;    
    }
    
    vector<vector<T>> res(r,vector<T>(c));
    
    for (int i = 0; i<r*c; i++){
        res[i/c][i%c] = mat[i/matCol][i%matCol];
    }
    
    return res;
}
template vector<vector<double>> matrix_tools<double>::
matrixReshape(vector<vector<double>> mat, int r, int c);
template vector<vector<float>> matrix_tools<float>::
matrixReshape(vector<vector<float>> mat, int r, int c);
template vector<vector<int>> matrix_tools<int>::
matrixReshape(vector<vector<int>> mat, int r, int c);

template <class T>
void matrix_tools<T>::save(vector<T> array , char fhead[], char path[], int n_realizations){

    FILE *output;

    char n1[] = "000", n2[] = "00", n3[] = "0";
    char ftail[]=".txt";
    char number1[15];
    char fullname[strlen(path)+strlen(ftail)+strlen(fhead)+strlen(number1)];

    mkdir(path);

    sprintf(number1,"%d", n_realizations);

    memset(fullname,'\0',strlen(path)+strlen(ftail)+strlen(fhead)+12+4);
    strcat(fullname,path);
    strcat(fullname,fhead);

    if (n_realizations<10) 
        strcat(fullname,n1);
    else if (n_realizations<100 && n_realizations>=10)
        strcat(fullname,n2);
    else if (n_realizations<1000 && n_realizations>=100)
        strcat(fullname,n3);
                        
    strcat(fullname,number1);
    strcat(fullname,ftail);

    output=fopen(fullname,"w");

    for (int i = 0; i < array.size(); i++){
        fprintf(output,"%d\t%.10lf\n", i, array[i]);
    }

    fclose(output);

}
template void matrix_tools<double>::save(vector<double> array , char fhead[], char path[], int n_realizations);
template void matrix_tools<float>::save(vector<float> array , char fhead[], char path[], int n_realizations);
template void matrix_tools<int>::save(vector<int> array , char fhead[], char path[], int n_realizations);
