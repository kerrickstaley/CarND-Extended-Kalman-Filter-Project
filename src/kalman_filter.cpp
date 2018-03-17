#include <cmath>
#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd y = z - H_ * x_;
  UpdateInner(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // assumes that the measurement function is the polar rho/gamma/rho_dot function used for radar
  VectorXd z_pred(3);
  double px = x_[0];
  double py = x_[1];
  double vx = x_[2];
  double vy = x_[3];

  if (px == 0 && py == 0) {
    z_pred << 0, 0, 0;
  } else {
    double rho = sqrt(px * px + py * py);

    z_pred <<
      rho,
      atan2(py, px),
      (vx * px + vy * py) / rho;
  }

  VectorXd y = z - z_pred;
  // normalize gamma component of y to be between -pi and pi
  // C++ modulo sign rule is that sign always matches the left operand, which makes this more complicated :(
  // TODO pull this out into a function and unit test it
  y[1] = fmod(fmod(y[1], 2 * M_PI) + 3 * M_PI, 2 * M_PI) - M_PI;

  UpdateInner(y);
}

void KalmanFilter::UpdateInner(const VectorXd &y) {
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();

  x_ = x_ + K * y;
  P_ -= K * H_ * P_;
}
