#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  if (!estimations.size()) {
    // no idea if throwing a string is a good practice, oh well
    throw std::string("empty vector passed to CalculateRMSE");
  }

  if (estimations.size() != ground_truth.size()) {
    throw std::string("`estimations` was a different size than `ground_truth`");
  }

  VectorXd sqDiffSum = VectorXd::Zero(estimations[0].size());

  for (int i = 0; i < estimations.size(); i++) {
    if (estimations[i].size() != sqDiffSum.size()) {
      throw std::string("inconsistently sized VectorXd passed in `estimations`");
    }
    if (ground_truth[i].size() != sqDiffSum.size()) {
      throw std::string("inconsistently sized VectorXd passed in `ground_truth`");
    }

    VectorXd diff = estimations[i] - ground_truth[i];
    diff = diff.array() * diff.array();
    sqDiffSum += diff;
  }

  VectorXd rmse = (sqDiffSum / estimations.size()).array().sqrt();
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  // calculates the Jacobian used in the Kalman filter measurement update for radar measurements
  // the Jacobian H is the partial derivative of each component of the measurement z = (rho, gamma, rho_dot) with
  // respect to the state x = (p_x, p_y, v_x, v_y), such that dz / dx = Hx (where z and x are column vectors)
  // Hence, H has 3 rows (one for each element of z) and 4 columns (one for each element of x)
  // Returns a zero matrix if division by zero occurs in calculation of result.
  MatrixXd H = MatrixXd::Zero(3, 4);

  double px = x_state[0];
  double py = x_state[1];
  double vx = x_state[2];
  double vy = x_state[3];

  double px2_py2 = px * px + py * py;
  double px2_py2_12 = sqrt(px2_py2);
  double px2_py2_32 = px2_py2 * px2_py2_12;

  // guard against division by zero
  if (px2_py2 == 0) {
    return H;
  }

  H <<
      // row 1, partial derivative of rho with respect to x
      // rho = sqrt(p_x^2 + p_y^2)
      px / px2_py2_12,
      py / px2_py2_12,
      0,
      0,
      // row 2, partial derivative of gamma with respect to x
      // gamma = atan(py / px)
      // note: d/dx atan(x) = 1 / (1 + x^2)
      - py / px2_py2,
      px / px2_py2,
      0,
      0,
      // row 3, partial derivative of rho_dot with respect to x
      // rho_dot = v . p / abs(v) = (vx * px + vy * py) / sqrt(px^2 + py^2)
      py * (vx * py - vy * px) / px2_py2_32,
      px * (vy * px - vx * py) / px2_py2_32,
      px / px2_py2_12,
      py / px2_py2_12;

  return H;
}
