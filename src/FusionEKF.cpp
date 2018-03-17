#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.Q_ = MatrixXd(4, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  // measurement function - laser
  H_laser_ <<
    1, 0, 0, 0,
    0, 1, 0, 0;

  // process covariance matrix
  // TODO: not sure what this should actually be
  ekf_.P_ = MatrixXd::Identity(4, 4);
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    ekf_.x_ = VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = measurement_pack.raw_measurements_[0];
      double gamma = measurement_pack.raw_measurements_[1];
      double rho_dot = measurement_pack.raw_measurements_[2];
      ekf_.x_ <<
        rho * cos(gamma),   // px
        rho * sin(gamma),   // py
        0,                  // vx
        0;                  // vy
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_ <<
        measurement_pack.raw_measurements_[0],  // px
        measurement_pack.raw_measurements_[1],  // py
        0,                                      // vx
        0;                                      // vy
    }

    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1e6;
  double dt2 = dt * dt;
  double dt3 = dt2 * dt;
  double dt4 = dt3 * dt;

  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.F_ <<
    1, 0, dt,  0,
    0, 1,  0, dt,
    0, 0,  1,  0,
    0, 0,  0,  1;

  int noise_ax = 9;
  int noise_ay = 9;

  ekf_.Q_ <<
      // row 1
      dt4 / 4.0 * noise_ax,
      0,
      dt3 / 2.0 * noise_ax,
      0,
      // row 2
      0,
      dt4 / 4.0 * noise_ay,
      0,
      dt3 / 2.0 * noise_ay,
      // row 3
      dt3 / 2.0 * noise_ax,
      0,
      dt2 * noise_ax,
      0,
      // row 4
      0,
      dt3 / 2.0 * noise_ay,
      0,
      dt2 * noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  } else {
    ekf_.R_ = R_radar_;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
