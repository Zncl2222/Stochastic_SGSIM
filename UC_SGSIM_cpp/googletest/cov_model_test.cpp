# include "../header/cov_model.hpp"

# include <gtest/gtest.h>
# include <vector>


using namespace std;

TEST(testCase, cov_model_getparameters1){

    Gaussian<double> covmodel(17.32, 1 ,35, 1);

    double range = covmodel.get_m_range();
    double sill = covmodel.get_m_sill();

    vector<double> res = {range, sill};
    vector<double> expect = {17.32, 1};
    EXPECT_EQ(res,expect);
}

TEST(testCase, cov_model_getparameters2){

    Gaussian<double> covmodel(17.32, 1 ,35, 1);


    int hs = covmodel.get_m_hs();
    int bw = covmodel.get_m_bw();

    vector<int> res = {hs, bw};
    vector<int> expect = {35, 1};

    EXPECT_EQ(res,expect);
}


int main(int argc, char** argv){

    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}