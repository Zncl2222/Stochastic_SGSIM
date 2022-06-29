# pragma once
# include <vector>
# include <cmath>

using namespace std;

template <class T>
class covariance_base{

    public:
        vector<T> vario_array;
        covariance_base(){};
        covariance_base(T range, T sill, int hs, int bw):
        m_range(range), m_sill(sill), m_hs(hs), m_bw(bw){
            this->vario_array = vector<T>(hs);
        };

        void variogram(vector<T> array, int mlen);

        vector<T> get_cov1d(){
            return this->cov1d;
        }
        vector<T> get_cov2d(){
            return this->cov2d;
        }

        void set_cov1d_size(int row_size){
            this->cov1d = vector<T>(row_size);
        }

        void set_cov2d_size(int row_size, int col_size){
            this->cov2d = vector<T>(row_size*col_size);
        }

        vector<T> get_vario(){
            return this->vario;
        }

        T get_m_range(){return this->m_range;}
        T get_m_sill(){return this->m_sill;}
        int get_m_hs(){return this->m_hs;}
        int get_m_bw(){return this->m_bw;}

    protected:
        vector<T> cov1d, cov2d;
        T m_range, m_sill;
        int m_hs;
        int m_bw;

}; 

template<class T>
class Gaussian: public covariance_base<T>{

    public:
        Gaussian():covariance_base<T>(){};
        Gaussian(T range, T sill, int hs, int bw):covariance_base<T>(range, sill, hs, bw){};

        void Gaussian_model(vector<T> arr, int n);

        void Gaussian_model(vector<vector<T>> arr, int n);
};