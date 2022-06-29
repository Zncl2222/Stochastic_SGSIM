# ifndef SGSIM_HPP
# define SGSIM_HPP

# include <cstdio>
# include <cmath>
# include <vector>
# include <iostream>
# include "../header/cov_model.hpp"

using namespace std;

template <class T> 
class sgsim_array{

    public:
        vector<vector<T>> randomfield;

        sgsim_array(int x_grid, int nR);
    
        void sgsim(Gaussian<T> covmodel, int randomseed = 9999);

        void print_randomfield(){
            
            for(int i = 0 ; i < this->randomfield.size(); i++){
                cout << "Realizations: " << i << endl;    
                for(int j = 0; j < this->randomfield[0].size(); j++){
                    cout << this->randomfield[i][j] << ", ";
                }
                cout << endl;
            }
        }

        int get_x_length(){return this->x_length;}
        int get_nR(){return this->nR;}

    private:
        int x_length, nR;


};

# endif