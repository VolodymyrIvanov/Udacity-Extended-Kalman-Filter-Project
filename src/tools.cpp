#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse = VectorXd::Zero(4);

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() == 0 || estimations.size() != ground_truth.size()) {
    return rmse;
  }

  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array().pow(2);
    rmse = rmse + residual;
  }

  //calculate the mean
  rmse = rmse.array() / estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj = MatrixXd::Zero(3,4);
  //recover state parameters
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);
  double pxy = px * px + py * py;

  //check division by zero
  if (pxy > THRESHOLD)
  {
    //compute the Jacobian matrix
    float h11 = px / sqrt(pxy);
    float h12 = py / sqrt(pxy);
    float h21 = -py / pxy;
    float h22 = px / pxy;
    float h31 = py * (vx * py - vy * px) / pow(pxy, 1.5);
    float h32 = px * (vy * px - vx * py) / pow(pxy, 1.5);
    float h33 = h11;
    float h34 = h12;
    Hj << h11, h12,   0,   0,
          h21, h22,   0,   0,
          h31, h32, h33, h34;
  }

  return Hj;
}

VectorXd Tools::TransformFromCartesianToPolar(const VectorXd& cartesian) {
  VectorXd polar = VectorXd(3);
  double px = cartesian[0];
  double py = cartesian[1];
  double vx = cartesian[2];
  double vy = cartesian[3];

  double rho = sqrt(px * px + py * py);
  double theta = atan2(py, px);
  double rhodot = 0.0;
  if (rho > THRESHOLD) {
    rhodot = (px * vx + py * vy) / rho;
  }

  polar << rho, theta, rhodot;
  return polar;
}

VectorXd Tools::TransformFromPolarToCartesian(const VectorXd& polar) {
  VectorXd cartesian = VectorXd(4);
  double rho = polar[0];
  double theta = polar[1];
  double rhodot = polar[2];

  double px = rho * cos(theta);
  double py = rho * sin(theta);
  double vx = rhodot * cos(theta);
  double vy = rhodot * sin(theta);

  cartesian << px, py, vx, vy;
  return cartesian;
}
