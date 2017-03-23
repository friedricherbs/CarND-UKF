#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

#include "measurement_package.h"
#include "ground_truth_package.h"

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
  * A helper method to estimate measurement noise.
  */
  void EstimateMeasurementNoise(const vector<MeasurementPackage>& measurements, const vector<GroundTruthPackage>& ground_truth);

  /**
  * A helper method to estimate process noise.
  */
  void EstimateProcessNoise(const vector<GroundTruthPackage>& ground_truth);
};

#endif /* TOOLS_H_ */
