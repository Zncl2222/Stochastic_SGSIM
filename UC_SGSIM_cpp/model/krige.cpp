# include <vector>
# include <cstdio>
# include <cmath>
# include "../header/krige.hpp"
# include "../header/random_tools.hpp"
# include "../header/sort_tools.hpp"
# include "../header/matrix_tools.hpp"
# include "../header/cov_model.hpp"

template <class T>
kriging<T>::kriging(Gaussian<T> covmodel){

    this->m_covmodel = covmodel;

    this->dist_temp = vector<T>(6);
    //this->flatten_temp = vector<T>(6*6);
    this->distcov_temp2 = vector<T>(6);
    //this->distcov_temp = vector<T>(6);
    this->weights = vector<T>(6);
    this->pdist_temp = vector<vector<T>>(6,vector<T>(6));
    this->datacov_temp = vector<vector<T>>(6,vector<T>(6));
    
}
template kriging<double>::kriging(Gaussian<double> covmodel);
template kriging<float>::kriging(Gaussian<float> covmodel);
template kriging<int>::kriging(Gaussian<int> covmodel);

template <class T>
void kriging<T>::simplekriging(vector<T>& array, vector<T> sampled, vector<T> u_array,
int currlen, T unsampled_point, int neighbor, int randomseed){

    random_tools<T> rt;
    T range = this->m_covmodel.get_m_range();
    T sill = this->m_covmodel.get_m_sill();
    this->array2d_temp = vector<vector<T>>(3,vector<T>(array.size()));

    if(neighbor == 0){
        
        array[(int)unsampled_point] = rt.random_normal(randomseed);
    }
    else{

        int close = 0;

        for(int j = 0; j < currlen; j++){
        
            u_array[j] = fabs(sampled[j] - unsampled_point);
            if(u_array[j] < range * 1.732){
                close++;
            }
        }

        if(close == 0){

            array[(int)unsampled_point] = rt.random_normal(randomseed);
        }
        else{

            for(int j = 0; j < currlen; j++){

                array2d_temp[0][j] = sampled[j];
                array2d_temp[1][j] = array[(int)sampled[j]];
                array2d_temp[2][j] = u_array[j];
            }

            if (neighbor >= 2){
                
                sorting<T> st;
                st.quicksorted_2d(array2d_temp, 0, currlen-1);
            }

            for(int j = 0; j < neighbor; j++){
                
                dist_temp[j] = array2d_temp[0][j];
                distcov_temp2[j] = array2d_temp[2][j]; 
            }
            
            matrix_tools<T> mt;
            
            mt.pdist(dist_temp, pdist_temp, neighbor);
            
            m_covmodel.Gaussian_model(pdist_temp, neighbor);
            
            flatten_temp = m_covmodel.get_cov2d();
            
            mt.matrixform(flatten_temp, datacov_temp, neighbor);
            
            m_covmodel.Gaussian_model(distcov_temp2, neighbor);

            distcov_temp = m_covmodel.get_cov1d();

            if(neighbor >= 1)
                mt.LUdecomposition(datacov_temp, distcov_temp, weights, neighbor); 

            T zvalue = 0;
            T svar = 0;
            T fix = 0;
            
            for(int j = 0; j < neighbor; j++){

                zvalue = zvalue + array2d_temp[1][j] * weights[j] * sill;
                svar = svar + distcov_temp[j] * weights[j];
            }

            svar = 1 - svar;
            
            if(svar < 0)
                svar = 0;
            fix = rt.random_normal(randomseed) * pow(svar, 0.5);
            
            array[(int)unsampled_point] = zvalue + fix;
        }
    }
}
template void kriging<double>::simplekriging(vector<double>& array, vector<double> sampled, vector<double> u_array,
int currlen, double unsampled_point, int neighbor, int randomseed);

template void kriging<float>::simplekriging(vector<float>& array, vector<float> sampled, vector<float> u_array,
int currlen, float unsampled_point, int neighbor, int randomseed);

template void kriging<int>::simplekriging(vector<int>& array, vector<int> sampled, vector<int> u_array,
int currlen, int unsampled_point, int neighbor, int randomseed);

