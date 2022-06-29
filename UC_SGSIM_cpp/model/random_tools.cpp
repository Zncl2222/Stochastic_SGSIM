# include "../header/random_tools.hpp"

using namespace std;
template<class T>
vector<T> random_tools<T>::randompath(vector<T> rpath, int length, int randomseed){ 

  srand(randomseed);
  int rtemp;
  for(int i=length-1;i;i--)
  {
    int rindex= rand()%i;
    
    do{
      rtemp=rpath[rindex];
      
      rpath[rindex]=rpath[i];
      
      rpath[i]=rtemp;

    }while(0);

  }

  return rpath;
}
template vector<double> random_tools<double>::randompath(vector<double> rpath, int length,int randomseed);
template vector<float> random_tools<float>::randompath(vector<float> rpath, int length,int randomseed);
template vector<int> random_tools<int>::randompath(vector<int> rpath, int length,int randomseed);


template<class T>
T random_tools<T>::random_normal(int randomseed){
    
    srand(randomseed);
    T x = rand() / (T)RAND_MAX;
    T y = rand() / (T)RAND_MAX;
    T z = sqrt(-2*log(x)) * cos(2*M_PI*y);

    return z;
}
template double random_tools<double>::random_normal(int randomseed);
template float random_tools<float>::random_normal(int randomseed);
template int random_tools<int>::random_normal(int randomseed);