#ifndef UKF_H
#define UKF_H
#include "Eigen/Dense"
#include "measurement_package.h"
#include "ground_truth_package.h"
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Initial state covariance terms
static const double V_INIT_VAR   = 0.2;
static const double PSI_INIT_VAR = 0.3;
static const double PSI_DOT_VAR  = 0.3F;

// Scaling parameters for Unscented transform
static const double ALPHA        = 0.1;
static const double BETA         = 2.0;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  // previous timestamp
  long long previous_timestamp_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  // Laser measurement covariance matrix
  MatrixXd R_laser_;

  // Radar measurement covariance matrix
  MatrixXd R_radar_;

  ///* Sigma points matrix
  MatrixXd Xsig_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   * @param gt_package The ground truth of the state x at measurement time
   */
  void ProcessMeasurement(const MeasurementPackage& meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(const double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(const MeasurementPackage& meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(const MeasurementPackage& meas_package);

private:

    /**
    * Initialize state from an initial radar measurement
    * @param meas_package The measurement at k+1
    */
    void InitRadar(const MeasurementPackage& meas_package);

    /**
    * Initialize state from an initial laser measurement
    * @param meas_package The measurement at k+1
    */
    void InitLaser(MeasurementPackage meas_package);

    /**
    * Create Sigma points for augmented state including process noise terms
    * @param Xsig_out Augmented Sigma points
    */
    void AugmentedSigmaPoints(MatrixXd* Xsig_out);

    /**
    * Predict Sigma points according to CTRV model
    * @param delta_t Time difference for prediction
    */
    void SigmaPointPrediction(const double delta_t);

    /**
    * Predict state mean and covariance matrix according to CTRV model
    */
    void PredictMeanAndCovariance();

    /**
    * Convert predicted sigma points to Radar measurement space, calculate mean predicted measurement and measurement covariance S
    * @param z_out    Predicted state in Radar measurement space
    * @param S_out    Measurement covariance matrix
    * @param Zsig_out Predicted sigma points in Radar measurement space
    */
    bool PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_out);

    /**
    * Convert predicted sigma points to Laser measurement space, calculate mean predicted measurement and measurement covariance S
    * @param z_out    Predicted state in Laser measurement space
    * @param S_out    Measurement covariance matrix
    * @param Zsig_out Predicted sigma points in Laser measurement space
    */
    void PredictLaserMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_out);
};

#endif /* UKF_H */
