
# include "../header/matrix_tools.hpp"


# include <gtest/gtest.h>
# include <vector>
# include <climits>

using namespace std;

TEST(testCase, arange_size_normal){
    matrix_tools<double> mt;
    vector<double> arr = mt.arange(150);

    EXPECT_EQ(arr.size(),150);
}

TEST(testCase, arange_size_zero){
    matrix_tools<double> mt;
    vector<double> arr = mt.arange(0);

    EXPECT_EQ(arr.size(),0);
}
/*
TEST(testCase, arange_size_INTMAX){
    
    vector<double> arr = arange<double>(INT_MAX);

    EXPECT_EQ(arr.size(),INT_MAX);
}*/

TEST(testCase, pdist_normal){
    matrix_tools<double> mt;
    vector<double> arr = {1,2,3};
    vector<vector<double>> expect = {{0,1,2},{1,0,1},{2,1,0}};
    vector<vector<double>> pdist_temp(3,vector<double>(3));
    mt.pdist(arr, pdist_temp,arr.size());
    EXPECT_EQ(pdist_temp, expect);
}

TEST(testCase, pdist_zero){

    matrix_tools<double> mt;
    vector<double> arr = {0};
    vector<vector<double>> expect = {{0}};
    vector<vector<double>> pdist_temp(1,vector<double>(1));
    mt.pdist(arr, pdist_temp, arr.size());
    EXPECT_EQ(pdist_temp, expect);
}


TEST(testCase, matrixform_normal){
    matrix_tools<double> mt;
    vector<double> arr = {0,1,2,3,4,5,6,7,8};
    vector<vector<double>> res(3,vector<double>(3));
    vector<vector<double>> expect = {{0,1,2},{3,4,5},{6,7,8}};
    mt.matrixform(arr, res, res.size());
    EXPECT_EQ(res, expect);
}

TEST(testCase, matrixform_none){

    matrix_tools<double> mt;
    vector<double> arr = {0};
    vector<vector<double>> res(1,vector<double>(1));
    vector<vector<double>> expect = {{0}};
    mt.matrixform(arr, res, res.size());
    EXPECT_EQ(res, expect);
}

TEST(testCase, matrixReshape_rowsize_normal){

    matrix_tools<double> mt;
    vector<vector<double>> arr(10,vector<double>(3));
    vector<vector<double>> res;
    res = mt.matrixReshape(arr, 5 , 6);
    EXPECT_EQ(res.size(),5);
}

TEST(testCase, matrixReshape_colsize_normal){
    matrix_tools<double> mt;
    vector<vector<double>> arr(10,vector<double>(3));
    vector<vector<double>> res;
    res = mt.matrixReshape(arr, 5 , 6);
    EXPECT_EQ(res[0].size(),6);
}

TEST(testCase, matrixReshape_normal){
    matrix_tools<double> mt;
    vector<vector<double>> arr = {{1,2,3},{4,5,6},{7,8,9},{10,11,12}};
    vector<vector<double>> expect = {{1,2,3,4,5,6},{7,8,9,10,11,12}};
    vector<vector<double>> res;
    res = mt.matrixReshape(arr, 2 , 6);
    EXPECT_EQ(res,expect);
}


int main(int argc, char** argv){

    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}