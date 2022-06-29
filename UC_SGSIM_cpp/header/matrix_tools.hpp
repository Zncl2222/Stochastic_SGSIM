# pragma once
# include <cstdio>
# include <cmath>
# include <vector>
# include <string.h>
//# include <sys/io.h> // for linux
//# include <sys/stat.h> // for linux
# include <io.h>     // for windows
using namespace std;

template<class T>
class matrix_tools{

    public:
        
        void LUdecomposition(vector<vector<T>> a, vector<T> b, vector<T>& x, int n);
        vector<T> arange(int x);
        void pdist(vector<T> x, vector<vector<T>>& c, int n_dim);
        void matrixform(vector<T> x,vector<vector<T>>& matrix, int n_dim);
        vector<vector<T>> matrixReshape(vector<vector<T>> mat, int r, int c);
        void save(vector<T> array , char fhead[], char path[], int n_realizations);

};
