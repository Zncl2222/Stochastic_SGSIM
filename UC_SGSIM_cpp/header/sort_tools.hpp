# ifndef SORT_TOOLS_HPP
# define SORT_TOOLS_HPP


# include <iostream>
# include <vector>
# include <cmath>


using namespace std;

template <class T> 
class sorting{

    public:

        static void swap2d(vector<vector<T>>& arr1, int i, int j);

        static int partition2d(vector<vector<T>>& matrix, int front, int end);

        static void quicksorted_2d(vector<vector<T>>& matrix, int front, int end);

};
# endif