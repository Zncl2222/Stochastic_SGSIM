# pragma once
# include <cstdio>
# include <vector>
# include <cmath>


using namespace std;

template <class T>
class random_tools{

    public:
        vector<T> randompath(vector<T> rpath, int length, int randomseed);
        T random_normal(int randomseed);

};
