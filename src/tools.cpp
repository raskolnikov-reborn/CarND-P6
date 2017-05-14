#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

/**
 * Constructor: does nothing
 */
Tools::Tools() {}

/**
 * Destructor: does nothing
 */
Tools::~Tools() {}

/**
 * CalculateRMSE: calculates the Root Mean Square Error to evaluate performance
 * @param estimations: filter outputs
 * @param ground_truth: ground truth values
 */
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
		const vector<VectorXd> &ground_truth)
{
	// Initialize Vector
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// Check for Valid inputs
	if (estimations.size() != ground_truth.size())
	{
		std::cout << " Size Mismatch!" <<std::endl;
		return rmse;
	}

	if ((estimations.size()==0) || (ground_truth.size()==0))
	{
		std::cout << "Empty Vector!" <<std::endl;
		return rmse;
	}

	// Accumulate square residual error
	for(int i=0; i < estimations.size(); ++i){
		// ... your code here

		VectorXd err = ground_truth[i] - estimations[i];
		VectorXd err_sq = err.array()*err.array();

		rmse += err_sq;

	}

	// Divide by size to get mean value
	rmse = rmse/estimations.size();

	// Find square root
	rmse = rmse.array().sqrt();


	return rmse;
}

/**
 * CalculateJacobian: calculates the Jacobian matrix of partial derivatives
 * @param x_state: State variable
 */
MatrixXd Tools::CalculateJacobian(const VectorXd& x_state)
{
	// Initialze the value
	MatrixXd Hj(3,4);

	//recover state parameters from argument
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	// Calculate Squared sum and its root
	float sum = px*px + py*py;
	float sum_sq = pow (sum,0.5);
	float sum_3_2 = sum*sum_sq;

	//check division by zero
	if (sum==0)
	{
		std::cout<<"Division by Zero error";
		return Hj;
	}

	// Compute the Elements for the Jacobian
	float ele_1_1 = px/sum_sq;
	float ele_1_2 = py/sum_sq;

	float ele_2_1 = -py/sum;
	float ele_2_2 = px/sum;

	float ele_3_1 = py*(vx*py - vy*px)/sum_3_2;
	float ele_3_2 = px*(vy*px - vx*py)/sum_3_2;

	float ele_3_3 = ele_1_1;
	float ele_3_4 = ele_1_2;


	// Update Jacobian matrix
	Hj <<   ele_1_1, ele_1_2,0,0,
			ele_2_1, ele_2_2,0,0,
			ele_3_1,ele_3_2,ele_3_3,ele_3_4;

	return Hj;
}


/**
 * getCovariance: calculates the Covariance Matrix
 * @param dt: timestep
 * @param noise_ax: Variance in x
 * @param noise_ay: Variance in y
 */
MatrixXd Tools::getCovariance(const double& dt, const float& noise_ax, const float& noise_ay)
{

	// Referring to sensor-fusion-ekf-reference.pdf
	// Section 11.4

	// Computing Q = GQvGt
	MatrixXd G = MatrixXd(4,2);

	float dt2 = dt*dt;


	G <<    dt2/2,0,
			0,dt2/2,
			dt,0,
			0,dt;

	MatrixXd Qv = MatrixXd(2,2);

	Qv <<   noise_ax,0,0,noise_ay;

	MatrixXd Q = G*Qv*G.transpose();

	return Q;
}


VectorXd Tools::cartesian2polar(const VectorXd& x_state) {
	VectorXd hx_prime(3);

	// Recover State Values
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	// Calculate radius
	float radius = sqrt(px*px + py*py);

	// Min cap to avoid divide by zero
	if (radius < 0.001) {
		radius = 0.001;
	}

	// Calculate Bearing
	float bearing = 0.0;
	// check for zero by zero and calculate tan inverse
	if (fpclassify(px) != FP_ZERO && fpclassify(py) != FP_ZERO) {
		bearing = atan2(py, px);
	}

	// Calculate radial velocity
	float delta_r = (px*vx + py*vy) / radius;

	// Convert to rho,theta,rho_dot
	hx_prime << radius, bearing, delta_r;

	return hx_prime;
}

VectorXd Tools::polar2cartesian(const Eigen::VectorXd & z)
{
	VectorXd result(4);

	// recover values
	float radius = z(0);
	float bearing = z(1);
	float delta_r = z(2);

	// calculate cosine and sine
	float cos_b = cos(bearing);
	float sin_b = sin(bearing);

	// Calculate Cartesian Space value
	float px = radius*cos_b;
	float py = radius*sin_b;
	float vx = delta_r*cos_b;
	float vy = delta_r*sin_b;

	result << px, py, vx, vy;

	return result;
}
