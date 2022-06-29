# include <cstdio>
# include <vector>
# include <iostream>
# include <ctime>
# include "../header/sgsim.hpp"
# include "../header/krige.hpp"
# include "../header/matrix_tools.hpp"
# include "../header/sort_tools.hpp"
# include "../header/cov_model.hpp"


using namespace std;


template<class T>
void print(vector<T> arr){

    for(auto i = arr.begin(); i != arr.end(); i++){

        printf("%.3f, ", *i);
    }
    cout << endl;
}

template<class T>
void print2d(vector<vector<T>>& a){

    for(auto itr_row = a.begin(); itr_row != a.end(); itr_row++){

        for(auto itr_col = itr_row->begin(); itr_col != itr_row->end(); itr_col++){

            printf("%.3f ", *itr_col);
        }
        cout << "\n";
    }

}

int main(){

    clock_t start = clock();

    Gaussian<double> gaussian(17.32, 1 ,35, 1);
    sgsim_array<double> x(150, 1000);

    x.sgsim(gaussian);
    clock_t end = clock();
    //x.print_randomfield();
    double duration = double(end-start)/CLK_TCK;

    cout << "Duration: " << duration <<"s"<<endl; 
    system("Pause");

    return 0;

}
