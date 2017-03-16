#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace Tools {
  /**
  * A helper method to calculate RMSE.
  */
  VectorXd CalculateRMSE(const std::vector<VectorXd> &estimations, const std::vector<VectorXd> &ground_truth);

  /**
  * A helper method to calculate Jacobians.
  */
  MatrixXd CalculateJacobian(const VectorXd& x_state);

  /**
  * A helper method to map the state vector for RADAR
  */
  VectorXd MapStateVectorToMeasurementSpace(const VectorXd& x_state);

};

#endif /* TOOLS_H_ */
