#include <gtest/gtest.h>
#include <cmath>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;


TEST(CalculateRMSETest, Case1) {
  vector<VectorXd> estimations;
  VectorXd e(2);
  e << 0.0, 1.0;
  estimations.push_back(e);
  e << 2.0, 3.0;
  estimations.push_back(e);
  e << 4.0, 5.0;
  estimations.push_back(e);

  vector<VectorXd> ground_truths;
  VectorXd gt(2);
  gt << 4.0, 2.0;
  ground_truths.push_back(gt);
  gt << 2.0, 4.0;
  ground_truths.push_back(gt);
  gt << 0.0, 6.0;
  ground_truths.push_back(gt);

  Tools t;
  VectorXd rmse = t.CalculateRMSE(estimations, ground_truths);

  ASSERT_EQ(rmse.size(), 2);
  EXPECT_FLOAT_EQ(rmse[0], sqrt(32.0 / 3));
  EXPECT_FLOAT_EQ(rmse[1], 1);
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
