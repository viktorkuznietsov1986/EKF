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

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  


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
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);

    previous_timestamp_ = measurement_pack.timestamp_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float ro = measurement_pack.raw_measurements_[0];
      float fi = measurement_pack.raw_measurements_[1];

      float px = ro*sin(fi);
      float py = ro*cos(fi);

      ekf_.x_ << px, py, 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

   //compute the time elapsed between the current and previous measurements
	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = measurement_pack.timestamp_;

  //1. Modify the F matrix so that the time is integrated
	ekf_.F_(0,2) = dt;
	ekf_.F_(1,3) = dt;

  //2. Set the process covariance matrix Q
  int noise_ax = 9;
  int noise_ay = 9;
  float dt2 = dt*dt;
  float dt3 = dt2*dt/2;
  float dt4 = dt3*dt/2;

	ekf_.Q_ = MatrixXd(4,4);
	ekf_.Q_ << dt4*noise_ax, 0, dt3*noise_ax, 0,
	          0, dt4*noise_ay, 0, dt3*noise_ay,
	          dt3*noise_ax, 0, dt2*noise_ax, 0,
	          0, dt3*noise_ay, 0, dt2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */


  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    VectorXd z(3);

    z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], measurement_pack.raw_measurements_[2];
    //ekf_.UpdateEKF(z);
  } else {
    // Laser updates
    VectorXd z(2);
	  z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1];
    ekf_.Update(z);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}