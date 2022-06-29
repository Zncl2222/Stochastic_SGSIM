# include <cstdio>
# include <cmath>
# include <vector>
# include "cov_model.hpp"
# include "sgsim.hpp"

using namespace std;

template <class T>
class kriging{

    public:
        Gaussian<T> m_covmodel;
        kriging(){};
        kriging(Gaussian<T> covmodel);

        void simplekriging(vector<T>& array, vector<T> sampled, vector<T> u_array, 
                            int currlen, T unsampled_point, int neighbor, int randomseed);


    private:
        vector<T> array, sampled, u_array;
        int currlen, neighbor, randomseed;
        T unsampled_point;
        
        vector<T> dist_temp, flatten_temp, distcov_temp2, distcov_temp, weights;
        vector<vector<T>> pdist_temp, datacov_temp, array2d_temp;

};