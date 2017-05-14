#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
		MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in)
{
	x_ = x_in;
	P_ = P_in;
	F_ = F_in;
	H_ = H_in;
	R_ = R_in;
	Q_ = Q_in;
}


/**
 * Predict: Prediction step of the Kalman Filter
 */
void KalmanFilter::Predict()
{

	// Referring to sensor-fusion-ekf-reference.pdf
	// Section 5.1

	//x′ =Fx+u (11)
	x_ = F_*x_;

	//P′=FPFT+Q (12)
	P_ = F_*P_*F_.transpose() + Q_;

}
/**
 * Update: Update step of the Kalman Filter
 * @param z: raw measurement value
 */
void KalmanFilter::Update(const VectorXd &z)
{

	// Referring to sensor-fusion-ekf-reference.pdf
	// Section 5.1

	//	y = z − Hx′ (13)
	VectorXd y = z - H_ * x_;

	//	S=HP′HT+R (14)
	MatrixXd Ht = H_.transpose();
	MatrixXd PHt = P_ * Ht;
	MatrixXd S = H_ * PHt + R_;

	//	K = P′HT S−1 (15)
	MatrixXd K = PHt * S.inverse();

	//	x = x′ + Ky (16)
	x_ = x_ + K*y;

	//	P=(I−KH)P′ (17)
	MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
	P_ = (I - K*H_)*P_;

}


/**
 * Update: Update step of the Extended Kalman Filter
 * @param z: raw measurement value
 */
void KalmanFilter::UpdateEKF(const VectorXd &z)
{
	// Instantiate Tools object
	Tools tools;

	// Convert cartesian update to polar coordinates
	MatrixXd hx = tools.cartesian2polar(x_);


	// Referring to sensor-fusion-ekf-reference.pdf
	// Section 5.1

	//	y = z − Hx′ (13)
	VectorXd y = z - hx;

	while (y[1] > M_PI)
	{
		y[1] -= 2*M_PI;
	}

	while (y[1] < -M_PI)
	{
		y[1] += 2*M_PI;
	}

	//	S=HP′HT+R (14)
	MatrixXd Ht = H_.transpose();
	MatrixXd PHt = P_ * Ht;
	MatrixXd S = H_ * PHt + R_;

	//	K = P′HT S−1 (15)
	MatrixXd K = PHt * S.inverse();

	//	x = x′ + Ky (16)
	x_ = x_ + K*y;

	//	P=(I−KH)P′ (17)
	MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
	P_ = (I - K*H_)*P_;
}
