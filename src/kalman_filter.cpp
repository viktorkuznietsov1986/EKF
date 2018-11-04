#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

#include <iostream>
using namespace std;

// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

const double epsilon = 0.0001;

const double pi = 3.14159265358979323846;

KalmanFilter::KalmanFilter() {
}

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
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd y = z - H_ * x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  // introduce nonlinearity
  VectorXd hx = GetHx();
  VectorXd y = z - hx;
  normalizeBearing(y);

  // approximate the nonlinearity by using the corresponding Jacobian matrix
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::normalizeBearing(VectorXd& y) {
  double fi = y(1);

  // normalize to have fi in -pi:pi.
  while (fi < -pi || fi > pi) {
    if (fi < -pi) {
      fi += 2.0*pi;
    }
    else {
      fi -= 2.0*pi;
    }
  }

  y(1) = fi;
}

VectorXd KalmanFilter::GetHx() {
  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];

  VectorXd hx = VectorXd(3);

  hx << 0.0, 0.0, 0.0;

  if (fabs(px) < epsilon) {
    cout << "KalmanFilter::GetHx - px == 0. Can't calculate Hx." << endl;
    return hx;
  }

  // get range
  float ro = sqrt(px*px + py*py);

  if (fabs(ro) < epsilon) {
    cout << "KalmanFilter::GetHx - ro == 0. Can't calculate Hx." << endl;
    return hx;
  }

  // get bearing
  float fi = atan2(py, px);

  // get radial velocity
  float ro_dot = (px*vx + py*vy) / ro;

  // set up the non-linear measurement update vector
  hx << ro, fi, ro_dot;

  return hx;
}
