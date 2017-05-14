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

  // Laser Measurement matrix (2x4)
  H_laser_ << 1, 0, 0, 0,
		  	  0, 1, 0, 0;

  // Radar Measurement matrix (4x4)
  Hj_ <<  0, 0, 0, 0,
		  0, 0, 0, 0,
		  0, 0, 0, 0;

  // Process noise
  noise_ax_ = 9.0;
  noise_ay_ = 9.0;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_)
  {
	  // Save the current timestamp into previous
	  previous_timestamp_ = measurement_pack.timestamp_;
	  // first measurement
	  cout << "EKF: " << endl;
	  ekf_.x_ = VectorXd(4);
	  ekf_.x_ << 1, 1, 1, 1;


	  // Initialize an empty Process Covariance matrix
	  ekf_.Q_ = MatrixXd(4, 4);

	  // Initialize the State state transition matrix with dt=1
	  ekf_.F_ = MatrixXd(4, 4);
	  ekf_.F_ <<  1, 0, 1, 0,
				  0, 1, 0, 1,
				  0, 0, 1, 0,
				  0, 0, 0, 1;

	  // Initialize an empty Prediction Covariance matrix
	  ekf_.P_ = MatrixXd(4, 4);


	  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
	  {
		  /**
      Convert radar from polar to cartesian coordinates and initialize state.
		   */

		  // Initial Covariance
		  // velocity variance should be lower than position
		  // NOTE: Tuning might effect RMS
		  ekf_.P_ <<  5, 0, 0, 0,
				  	  0, 5, 0, 0,
				  	  0, 0, 0.2, 0,
				  	  0, 0, 0, 0.2;


		  // Convert Radar Polar Readings to Cartesian Values
		  ekf_.x_ = tools.polar2cartesian(measurement_pack.raw_measurements_);
	  }
	  else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
	  {
		  /**
      Initialize state.
		   */
		  // Initial Covariance
		  // velocity variance should be higher than position
		  // NOTE: Tuning might effect RMS
		  ekf_.P_ <<  0.2, 0, 0, 0,
				  	  0, 0.2, 0, 0,
				  	  0, 0, 10, 0,
				  	  0, 0, 0, 10;

		  // Extract laser data
		  double x = measurement_pack.raw_measurements_[0];
		  double y = measurement_pack.raw_measurements_[1];
		  ekf_.x_ <<  x,
				  	  y,
				  	  0.0,
				  	  0.0;
	  }

	  // done initializing, no need to predict or update
	  is_initialized_ = true;
	  return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  // Calculate elapsed time
  // Timestamp is in microseconds
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  // High Pass filter for timestep
  if (dt > 0.01) {
	  // Update State Transition Matrix
	  ekf_.F_(0, 2) = dt;
	  ekf_.F_(1, 3) = dt;

	  // Update Q
	  ekf_.Q_ = tools.getCovariance(dt, noise_ax_, noise_ay_);

	  ekf_.Predict();
  }

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
	  // Radar updates
	  // Calculate Jacobian
	  ekf_.H_ = tools.CalculateJacobian(ekf_.x_);

	  // Update Measurement update matrix
	  ekf_.R_ = R_radar_;

	  // EKF update
	  ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  }
  else
  {
	  // Laser updates

	  // Set measurement matrix
	  ekf_.H_ = H_laser_;

	  // Measurement update
	  ekf_.R_ = R_laser_;

	  // Run normal KF update
	  ekf_.Update(measurement_pack.raw_measurements_);

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
