#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

class Tools {
public:
	/**
	 * Constructor.
	 */
	Tools();

	/**
	 * Destructor.
	 */
	virtual ~Tools();

	/**
	 * A helper method to calculate RMSE.
	 */
	Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);

	/**
	 * A helper method to calculate Jacobians.
	 */
	Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd& x_state);

	/**
	 * A helper method to convert Polar to Cartesian Co-ordinates
	 */
	Eigen::VectorXd cartesian2polar(const Eigen::VectorXd& x_state);

	/**
	 * A helper method to convert Cartesian to Polar Co-ordinates
	 */
	Eigen::VectorXd polar2cartesian(const Eigen::VectorXd & z);


	/**
	 * A helper method to calculate the covariance matrix
	 */
	Eigen::MatrixXd getCovariance(const double& dt, const float& noise_ax, const float& noise_ay);




};

#endif /* TOOLS_H_ */
