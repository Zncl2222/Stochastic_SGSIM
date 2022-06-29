# include "../header/sgsim.hpp"
# include "../header/random_tools.hpp"
# include "../header/krige.hpp"
# include "../header/matrix_tools.hpp"
# include "../header/cov_model.hpp"
# include <iostream>

using namespace std;

template <class T>
sgsim_array<T>::sgsim_array(int x_length, int nR){
    this->x_length = x_length;
    this->nR = nR;
    this->randomfield = vector<vector<T>>(nR, vector<T>(x_length));
}
template sgsim_array<double>::sgsim_array(int x_length, int nR);
template sgsim_array<float>::sgsim_array(int x_length, int nR);
template sgsim_array<int>::sgsim_array(int x_length, int nR);

template <class T>
void sgsim_array<T>::sgsim(Gaussian<T> covmodel, int randomseed){

    kriging<T> krg(covmodel);
    matrix_tools<T> mt;
    random_tools<T> rt;

    vector<T> unsampled_array(this->x_length);
    vector<T> model_x(this->x_length);
    vector<T> x_randompath;
    vector<T> sampled_array(this->x_length);
    
    int count = 0;
    int currentlen;
    int neighbor;
    int flag;

    x_randompath = mt.arange(this->x_length);

    while(count < this->nR){

        cout<<"Realizations"<< count <<endl;
        currentlen = 0;
        neighbor = 0;
        flag = 0;

        x_randompath = rt.randompath(x_randompath, this->x_length, randomseed);

        for(int i = 0; i < this->x_length; i++){
            model_x[i] = 0;
            sampled_array[i] = 0;
            unsampled_array[i] = -1;
        }

        for(int i = 0; i < this->x_length; i++){

            krg.simplekriging(model_x, sampled_array, unsampled_array, currentlen, x_randompath[i], neighbor, randomseed);
            
            if(neighbor < 5)
                ++neighbor;
            
            sampled_array[i] = x_randompath[i];
            
            ++currentlen;
            ++randomseed;

            if (isfinite(model_x[x_randompath[i]]) == 0)
                ++flag;
            else
                this->randomfield[count][x_randompath[i]] = model_x[x_randompath[i]];
            
        }
        ++count;
        ++randomseed;

        if(flag == 0){
            char r_name[] = "Realizations";
            char r_path[] = "./Realizations/";
            char v_name[] = "Variogram";
            char v_path[] = "./Variogram/";
            mt.save(model_x, r_name, r_path ,count);

            covmodel.variogram(model_x, this->x_length);
            mt.save(covmodel.vario_array, v_name, v_path, count);
        }
        else
            count--;
    }
}
template void sgsim_array<double>::sgsim(Gaussian<double> covmodel, int randomseed);
template void sgsim_array<float>::sgsim(Gaussian<float> covmodel, int randomseed);
template void sgsim_array<int>::sgsim(Gaussian<int> covmodel, int randomseed);