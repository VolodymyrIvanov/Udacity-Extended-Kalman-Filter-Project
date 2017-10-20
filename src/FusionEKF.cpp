#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

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

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  //measurement matrix
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  //state covariance matrix P
  P_ = MatrixXd(4, 4);
  P_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1000, 0,
        0, 0, 0, 1000;

  //the initial transition matrix F_
  F_ = MatrixXd(4, 4);
  F_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;

  //the process covariance matrix
  Q_ = MatrixXd::Zero(4, 4);

  noise_ax = 9;
  noise_ay = 9;
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
    // first measurement
    cout << "EKF: " << endl;
    VectorXd x_ = VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      x_ << tools.TransformFromPolarToCartesian(measurement_pack.raw_measurements_);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      float px = measurement_pack.raw_measurements_[0];
      float py = measurement_pack.raw_measurements_[1];
      x_ << px, py, 0, 0;
    }

    ekf_.Init(x_, P_, F_, Q_);

    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  //Modify matrices
  F_(0, 2) = dt;
  F_(1, 3) = dt;

  Q_ << pow(dt, 4) * noise_ax / 4.0, 0, pow(dt, 3) * noise_ax / 2.0, 0,
        0, pow(dt, 4) * noise_ay / 4.0, 0, pow(dt, 3) * noise_ay / 2.0,
        pow(dt, 3) * noise_ax / 2.0, 0, pow(dt, 2) * noise_ax, 0,
        0, pow(dt, 3) * noise_ay / 2.0, 0, pow(dt, 2) * noise_ay;

  //Set modifications in Kalman Filter
  ekf_.BeforePredict(F_, Q_);
  //Prediction call
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    //Transform polar coordinates to Cartesian
    VectorXd state = tools.TransformFromPolarToCartesian(measurement_pack.raw_measurements_);
    //Evaluate Jacobian
    Hj_ = tools.CalculateJacobian(state);
    //Transform previous state from Cartesian to polar
    VectorXd Hx = tools.TransformFromCartesianToPolar(ekf_.x_);
    //Update call
    ekf_.UpdateEKF(measurement_pack.raw_measurements_, Hx, Hj_, R_radar_);
  } else {
    // Laser updates
    ekf_.Update(measurement_pack.raw_measurements_, H_laser_, R_laser_);
  }
}
