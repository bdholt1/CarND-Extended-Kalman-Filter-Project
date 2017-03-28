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
  // Hand-tuned the R_laser_ matrix to find good results.
  R_laser_ << 0.0225, 0,
              0, 0.0225;
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  R_radar_ = MatrixXd(3, 3);
  // Hand-tuned the R_radar_ matrix to find good results.
  R_radar_ << 0.9, 0, 0,
              0, 0.00009, 0,
              0, 0, 0.09;

  //the initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
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
    cout << "Kalman Filter Initialization " << endl;
    previous_timestamp_ = measurement_pack.timestamp_;

    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.P_ = MatrixXd(4, 4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      cout << "Initial measurement is RADAR" << endl;
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float rate = measurement_pack.raw_measurements_[2];
      ekf_.x_ << rho*cos(phi), rho*sin(phi), rate*cos(phi), rate*sin(phi);
      // If the initial measurement is radar then we are more confident of the velocity
      // and less confident of the position.
      ekf_.P_ << 1000, 0, 0, 0,
                 0, 1000, 0, 0,
                 0, 0, 1, 0,
                 0, 0, 0, 1;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      cout << "Initial measurement is LASER" << endl;
      /**
      Initialize state.
      */
      float px = measurement_pack.raw_measurements_[0];
      float py = measurement_pack.raw_measurements_[1];
      ekf_.x_ << px, py, 0, 0;
      // If the initial measurement is laser then we are more confident of the position
      // and have no confidence (therefore high uncertainty) in the velocity.
      ekf_.P_ << 1, 0, 0, 0,
                 0, 1, 0, 0,
                 0, 0, 1000, 0,
                 0, 0, 0, 1000;
    }

    cout << "initial state " << endl;
    cout << ekf_.x_ << endl;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
  //cout << "Predicting position to next timestamp" << endl;
  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;
  //cout << ekf_.F_ << endl;


  float dt4 = dt*dt*dt*dt;
  float dt3 = dt*dt*dt;
  float dt2 = dt*dt;


  float noise_ax = 4;
  float noise_ay = 4;
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << (dt4/4*noise_ax),                0,    (dt3/2*noise_ax),               0,
                            0,   (dt4/4*noise_ay),                 0,  (dt3/2*noise_ay),
             (dt3/2*noise_ax),                0,      (dt2*noise_ax),               0,
                            0,   (dt3/2*noise_ay),                 0,    (dt2*noise_ay);

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    cout << "Measurement is RADAR" << endl;
    // Radar updates
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    cout << "Measurement is LASER" << endl;
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
