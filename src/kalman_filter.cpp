#include "kalman_filter.h"


using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
    x_ = F_ * x_;
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z, const MatrixXd &H, const VectorXd &Hx, const MatrixXd &R) {
  
    VectorXd y = z - Hx;
    MatrixXd Ht = H.transpose();
    MatrixXd S = H * P_ * Ht + R;
    MatrixXd Si = S.inverse();
    MatrixXd K =  P_ * Ht * Si;
    x_ = x_ + (K * y);
    
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H) * P_;
}

//void KalmanFilter::UpdateRadar(const VectorXd &z, const MatrixXd &H, const MatrixXd &Hx, const MatrixXd &R) {
// 
//  VectorXd y = z - Hx;
//  MatrixXd Ht = H.transpose();
//  MatrixXd S = H * P_ * Ht + R;
//  MatrixXd Si = S.inverse();
//  MatrixXd K =  P_ * Ht * Si;
//  x_ = x_ + (K * y);
//  
//  long x_size = x_.size();
//  MatrixXd I = MatrixXd::Identity(x_size, x_size);
//  P_ = (I - K * H) * P_;
//
//  
//}
