# include "../header/sort_tools.hpp"

# include <gtest/gtest.h>
# include <vector>
# include <climits>

using namespace std;

TEST(testCase, quicksort2d_normal){

    sorting<double> st;

    vector<vector<double>> arr = {{1,2,9},
                                  {6,5,4},
                                  {9,8,2}};
    st.quicksorted_2d(arr,0,2);
    vector<vector<double>> expect = {{9,2,1},
                                     {4,5,6},
                                     {2,8,9}};
    EXPECT_EQ(arr,expect);
}

TEST(testCase, quicksort2d_zero){

    sorting<double> st;
    
    vector<vector<double>> arr = {{0},
                                  {0},
                                  {0}};
    st.quicksorted_2d(arr,0,2);
    vector<vector<double>> expect = {{0},
                                     {0},
                                     {0}};
    EXPECT_EQ(arr,expect);

}

TEST(testCase, quicksort2d_noraml_nochange){

    sorting<double> st;
    
    vector<vector<double>> arr = {{0,2147483647,151},
                                  {0,7,9},
                                  {0,1,2}};
    st.quicksorted_2d(arr,0,2);
    vector<vector<double>> expect = {{0,2147483647,151},
                                     {0,7,9},
                                     {0,1,2}};
    EXPECT_EQ(arr,expect);

}

int main(int argc, char** argv){

    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}