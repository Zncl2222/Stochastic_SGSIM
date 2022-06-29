# include "../header/cov_model.hpp"
# include "../header/matrix_tools.hpp"
# include <cmath>
# include <iostream>
# include <cstdio>

using namespace std;

template <class T>
void covariance_base<T>::variogram(vector<T> array, int mlen){

    T Z_temp;
    int count;
    vector<T> temp;
    vector<vector<T>> pdist_temp(mlen,vector<T>(mlen));
    matrix_tools<T> mt;
    
    temp = mt.arange(mlen);

    mt.pdist(temp,pdist_temp, mlen);

    for (int i = 0; i < this->m_hs; i += this->m_bw){   

        Z_temp = 0;
        count = 0;

        for (int j = 0; j < mlen; j++){

            for (int k = j + 1; k < mlen; k++){   
                
                if (pdist_temp[j][k] >= i - this->m_bw && pdist_temp[j][k] <= i + this->m_bw){

                    Z_temp = Z_temp + pow((array[j]-array[k]),2);
                    count += 1;
                    
                }
            }
        }
        
        if (Z_temp >= 1e-6)
            this->vario_array[i] = Z_temp / (2*count);
    }

}
template void covariance_base<double>::variogram(vector<double> mat, int mlen);
template void covariance_base<float>::variogram(vector<float> mat, int mlen);
template void covariance_base<int>::variogram(vector<int> mat, int mlen);

template <class T>
void Gaussian<T>::Gaussian_model(vector<T> arr, int n){
    
    this->set_cov1d_size(n);

    for(int i=0;i<n;i++){   

        this->cov1d[i]=this->m_sill-(this->m_sill*(1-exp(-3*pow(arr[i],2)/pow(this->m_range,2)))); 
    }
}
template void Gaussian<double>::Gaussian_model(vector<double> arr, int n);
template void Gaussian<float>::Gaussian_model(vector<float> arr, int n);
template void Gaussian<int>::Gaussian_model(vector<int> arr, int n);

template <class T>
void Gaussian<T>::Gaussian_model(vector<vector<T>> arr, int n){

    this->set_cov2d_size(n, n);

    for(int i=0;i<n;i++){

        for (int j=0;j<n;j++){

            this->cov2d[n*i+j]=this->m_sill-(this->m_sill*(1-exp(-3*pow(arr[i][j],2)/pow(this->m_range,2))));
        }
    }
}
template void Gaussian<double>::Gaussian_model(vector<vector<double>> arr, int n);
template void Gaussian<float>::Gaussian_model(vector<vector<float>> arr, int n);
template void Gaussian<int>::Gaussian_model(vector<vector<int>> arr, int n);
