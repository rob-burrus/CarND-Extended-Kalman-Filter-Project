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

  //measurement covariance matrix - laser
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  //measurement matrix
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  // State covariance
  P_ = MatrixXd(4, 4);
  P_ << 1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1000, 0,
      0, 0, 0, 1000;

  F_ = MatrixXd::Identity(4, 4);
  
  //Process Covariance
  Q_ = MatrixXd::Zero(4, 4);

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
   
    VectorXd x_in = VectorXd(4);
    previous_timestamp_ = measurement_pack.timestamp_;
    

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
      double rho_dot = measurement_pack.raw_measurements_[2];
      double x = rho * cos(phi);
      double y = rho * sin(phi);
      double x_dot = rho_dot * cos(phi);
      double y_dot = rho_dot * sin(phi);
      x_in << x, y, x_dot, y_dot;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      
      x_in << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }
    
    ekf_.Init(x_in, P_, F_, Q_);
    
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  //Update the state transition matrix F according to the new elapsed time.
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  
  float noise_ax = 9;
  float noise_ay = 9;
  
  //set the process covariance matrix Q
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
              0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
              dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
              0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
  
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
//    cout << "Radar " << endl;
    MatrixXd Hj_ = tools.CalculateJacobian(ekf_.x_);
    VectorXd Hx = tools.CartesianToPolar(ekf_.x_);
    ekf_.Update(measurement_pack.raw_measurements_, Hj_, Hx, R_radar_);
//    cout << "Radar 2" << endl;
  } else {
    // Laser updates
//    cout << "Laser " << endl;
    VectorXd Hx = H_laser_ * ekf_.x_;
    ekf_.Update(measurement_pack.raw_measurements_, H_laser_, Hx, R_laser_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
