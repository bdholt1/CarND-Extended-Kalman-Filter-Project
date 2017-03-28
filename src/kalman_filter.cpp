#include "kalman_filter.h"
#include "tools.h"
#include <iostream>

using namespace std;

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
  x_ = F_ * x_;  // TODO: Add support for B and u (control matrix and vector)
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd y = z - H_ * x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K =  P_ * Ht * S.inverse();

  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  x_ = x_ + (K * y);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  cout << "z = " << z << endl;
  VectorXd h_x = Tools::MapStateVectorToMeasurementSpace(x_);
  cout << "h_x = " << h_x << endl;
  VectorXd y = z - h_x;
  cout << "y = " << y << endl;
  if (y(1) >  M_PI)
    y(1) -= M_PI;
  else if (y(1) <  -M_PI)
    y(1) += M_PI;

  MatrixXd H_ = Tools::CalculateJacobian(x_);
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K =  P_ * Ht * S.inverse();

  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  x_ = x_ + (K * y);
  P_ = (I - K * H_) * P_;
}
