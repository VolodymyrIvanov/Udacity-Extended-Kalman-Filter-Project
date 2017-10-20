#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

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
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  /**
  * A helper method to calculate Jacobians.
  */
  MatrixXd CalculateJacobian(const VectorXd& x_state);

  /**
  * A helper method to convert coordinates from Cartesian to Polar.
  */
  VectorXd TransformFromCartesianToPolar(const VectorXd& cartesian);

  /**
  * A helper method to convert coordinates from Polar to Cartesian.
  */
  VectorXd TransformFromPolarToCartesian(const VectorXd& polar);
private:
  const double THRESHOLD = 0.0001;
};

#endif /* TOOLS_H_ */
