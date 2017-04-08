#include <iostream>
#include "ukf.h"

/**
 * Initializes Unscented Kalman filter
 */

#define CHECK_INNOVATION 0 // Skip update if innovation is too large
#define SUKF 0 // Use scaled Unscented transform

UKF::UKF() {

  // state vector and covariance matrix have not been set already
  is_initialized_ = false;

  //set state dimension
  n_x_ = 5;

  //set augmented dimension
  n_aug_ = 7;

  // init previous timestamp
  previous_timestamp_ = 0;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(n_x_);
  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_ << 100000, 0,     0,     0,      0,
        0,      100000,0,     0,      0,
        0,      0,     100000,0,      0,
        0,      0,     0,     100000, 0,
        0,      0,     0,     0,      100000; 

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5; 

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5; 

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15; 

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15; 

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3; 

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03; 

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3; 

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_*std_laspx_, 0,
              0,                     std_laspy_*std_laspy_;

  // initializing matrices
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_*std_radr_, 0,                       0,
              0,                   std_radphi_*std_radphi_, 0,
              0,                   0,                       std_radrd_*std_radrd_;

  // Sigma points matrix
  Xsig_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Weight vector for sigma points
  weights_ = VectorXd(2*n_aug_+1);
  weights_.fill(0.0);

  //define spreading parameter
  lambda_ = 3 - n_aug_;

  //current NIS for radar
  NIS_radar_ = 0.0;

  //current NIS for radar
  NIS_laser_ = 0.0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(const MeasurementPackage& meas_package) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    if ((meas_package.sensor_type_ == MeasurementPackage::RADAR) && use_radar_) 
    {
      // Init state from radar measurement
      InitRadar(meas_package);
    }
	else if ((meas_package.sensor_type_ == MeasurementPackage::LASER) && use_laser_) 
    {
        // Init state from radar measurement
        InitLaser(meas_package);
	}

	previous_timestamp_ = meas_package.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  //Prediction step
  //compute the time elapsed between the current and previous measurements
  const double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;
  Prediction(dt);

  //Update step
  if ((meas_package.sensor_type_ == MeasurementPackage::RADAR) && use_radar_) {
    // Radar updates
    UpdateRadar(meas_package);
  } else if ((meas_package.sensor_type_ == MeasurementPackage::LASER) && use_laser_) {
    // Laser updates
    UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(const double delta_t) {
  SigmaPointPrediction(delta_t);
  PredictMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const MeasurementPackage& meas_package) {
    //set measurement dimension, lidar can measure x and y
    const int n_z = 2;

    //set vector for weights
    VectorXd weights = VectorXd(2*n_aug_+1);
#if SUKF
    const double weight_0 = lambda_/(ALPHA*ALPHA*(lambda_+n_aug_)) + 1.0/(ALPHA*ALPHA) - 1.0;
#else
    const double weight_0 = lambda_/(lambda_+n_aug_);
#endif
    weights(0) = weight_0;
    for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
#if SUKF
        const double weight = 0.5/(ALPHA*ALPHA*(lambda_+n_aug_));
#else
        const double weight = 0.5/(n_aug_+lambda_);
#endif
        weights(i) = weight;
    }

    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    //create matrix for predicted measurement covariance
    MatrixXd S = MatrixXd(n_z,n_z);
    //create vector for mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    PredictLaserMeasurement(&z_pred, &S, &Zsig);

    //create vector for incoming radar measurement
    VectorXd z = VectorXd(n_z);
    z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points

        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        // state difference
        VectorXd x_diff = Xsig_.col(i) - x_;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        Tc = Tc + weights(i) * x_diff * z_diff.transpose();
    }

    //Kalman gain K;
    MatrixXd S_inv   = MatrixXd(n_z,n_z);
    Eigen::FullPivLU<MatrixXd> S_invertible(S);
    assert(S_invertible.isInvertible());
    S_inv            = S_invertible.inverse();
    const MatrixXd K = Tc * S_inv;

    //residual
    VectorXd z_diff = z - z_pred;

    //NIS 
    NIS_laser_ = z_diff.transpose() * S_inv * z_diff;

#if CHECK_INNOVATION
    // Robust filtering: innovation check
    if (NIS_laser_ > 13.82)
       return;
#endif

    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const MeasurementPackage& meas_package) 
{
  //set measurement dimension, radar can measure r, phi, and r_dot
  const int n_z = 3;

  //set vector for weights
  VectorXd weights = VectorXd(2*n_aug_+1);
#if SUKF
  const double weight_0 = lambda_/(ALPHA*ALPHA*(lambda_+n_aug_)) + 1.0/(ALPHA*ALPHA) - 1.0;
#else
  const double weight_0 = lambda_/(lambda_+n_aug_);
#endif

  weights(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
#if SUKF
    const double weight = 0.5/(ALPHA*ALPHA*(lambda_+n_aug_));
#else
    const double weight = 0.5/(n_aug_+lambda_);
#endif
    weights(i) = weight;
  }

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  //create matrix for predicted measurement covariance
  MatrixXd S = MatrixXd(n_z,n_z);
  //create vector for mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  if (!PredictRadarMeasurement(&z_pred, &S, &Zsig))
      return;

  //create vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd S_inv = MatrixXd(n_z,n_z);
  Eigen::FullPivLU<MatrixXd> S_invertible(S);
  assert(S_invertible.isInvertible());
  S_inv          = S_invertible.inverse();
  const MatrixXd K = Tc * S_inv;

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //NIS 
  NIS_radar_ = z_diff.transpose() * S_inv * z_diff;

#if CHECK_INNOVATION
  // Robust filtering: innovation check
  if (NIS_radar_ > 16.27)
    return;
#endif

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
}

void UKF::InitRadar(const MeasurementPackage& meas_package)
{
    const double cos_theta = cos(meas_package.raw_measurements_[1]);
    const double sin_theta = sin(meas_package.raw_measurements_[1]);
    const double x_init    = meas_package.raw_measurements_[0] * cos_theta;
    const double y_init    = meas_package.raw_measurements_[0] * sin_theta;
    const double v_init    = std::fabs(meas_package.raw_measurements_[2]);

    const double vx_init   = meas_package.raw_measurements_[2] * cos_theta;
    const double vy_init   = meas_package.raw_measurements_[2] * sin_theta;
    const double psi_init  = atan2(vy_init, vx_init);

    // state vector contains px, py, v, psi, psi_dot
    x_ = VectorXd(n_x_);
    x_ << x_init, y_init, v_init, psi_init, 0.0; 

    // Calculate initial covariance matrix via error propagation

    // Derivative of state coordinates with respect to ro
    const double dxdro  = cos_theta;
    const double dydro  = sin_theta;

    // Derivative of state coordinates with respect to theta
    const double dxdtheta  = -meas_package.raw_measurements_[0] * sin_theta;
    const double dydtheta  = meas_package.raw_measurements_[0]  * cos_theta;

    // Derivative of state coordinates with respect to ro dot
    const double dxdrodot  = 0;
    const double dydrodot  = 0;

    // Do error propagation
    const double var_x  = dxdro*dxdro*R_radar_(0,0)   + dxdtheta*dxdtheta*R_radar_(1,1)   + dxdrodot*dxdrodot*R_radar_(2,2);
    const double var_y  = dydro*dydro*R_radar_(0,0)   + dydtheta*dydtheta*R_radar_(1,1)   + dydrodot*dydrodot*R_radar_(2,2);

    P_ << var_x, 0,     0,          0,            0,
          0,     var_y, 0,          0,            0,
          0,     0,     V_INIT_VAR, 0,            0,
          0,     0,     0,          PSI_INIT_VAR, 0,
          0,     0,     0,          0,            PSI_DOT_VAR;
}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out)
{
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  const Eigen::LLT<MatrixXd> lltCalc(P_aug);
  assert(lltCalc.info() != Eigen::NumericalIssue);
  const MatrixXd L = lltCalc.matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
#if SUKF
    Xsig_aug.col(i+1)        = x_aug + ALPHA*sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - ALPHA*sqrt(lambda_+n_aug_) * L.col(i);
#else
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
#endif
  }

  //write result
  *Xsig_out = Xsig_aug;
}

void UKF::SigmaPointPrediction(const double delta_t)
{
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  AugmentedSigmaPoints(&Xsig_aug);

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    const double p_x      = Xsig_aug(0,i);
    const double p_y      = Xsig_aug(1,i);
    const double v        = Xsig_aug(2,i);
    const double yaw      = Xsig_aug(3,i);
    const double yawd     = Xsig_aug(4,i);
    const double nu_a     = Xsig_aug(5,i);
    const double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p    = v;
    double yaw_p  = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }

  //write result
  Xsig_ = Xsig_pred;
}

void UKF::PredictMeanAndCovariance()
{
  //create vector for weights
  VectorXd weights = VectorXd(2*n_aug_+1);

  //create vector for predicted state
  VectorXd x_pred = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P_pred = MatrixXd(n_x_, n_x_);

  // set weights
#if SUKF
  const double weight_0 = lambda_/(ALPHA*ALPHA*(lambda_+n_aug_)) + 1.0/(ALPHA*ALPHA) - 1.0;
#else
  const double weight_0 = lambda_/(lambda_+n_aug_);
#endif
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
#if SUKF
      const double weight = 0.5/(ALPHA*ALPHA*(lambda_+n_aug_));
#else
      const double weight = 0.5/(n_aug_+lambda_);
#endif
    weights(i) = weight;
  }

  //predicted state mean
  x_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_pred = x_pred + weights(i) * Xsig_.col(i);
  }

  //predicted state covariance matrix
  P_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_.col(i) - x_pred;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_pred = P_pred + weights(i) * x_diff * x_diff.transpose() ;
  }

  //write result
  x_ = x_pred;
  P_ = P_pred;
}

bool UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_out)
{
  //set measurement dimension, radar can measure r, phi, and r_dot
  const int n_z = 3;

  //set vector for weights
  VectorXd weights = VectorXd(2*n_aug_+1);
#if SUKF
  const double weight_0 = lambda_/(ALPHA*ALPHA*(lambda_+n_aug_)) + 1.0/(ALPHA*ALPHA) - 1.0;
#else
  const double weight_0 = lambda_/(lambda_+n_aug_);
#endif
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  
#if SUKF
      const double weight = 0.5/(ALPHA*ALPHA*(lambda_+n_aug_));
#else
      const double weight = 0.5/(n_aug_+lambda_);
#endif
    weights(i) = weight;
  }

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points

    // extract values for better readability
    const double p_x = Xsig_(0,i);
    const double p_y = Xsig_(1,i);
    const double v   = Xsig_(2,i);
    const double yaw = Xsig_(3,i);

    const double v1 = cos(yaw)*v;
    const double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);          //r
    Zsig(1,i) = atan2(p_y,p_x);                   //phi

    if (std::fabs(Zsig(0,i)) < 0.0001)
    {
        // Workaround: limit of px and py seems to be finite, 
        // probably close to zero as long as v1 and v2 are close to zero
        // TODO: use better approximation
        Zsig(2,i) = 0.0; 
        return false;
    }
    else
    {
        Zsig(2,i) = (p_x*v1 + p_y*v2 ) / Zsig(0,i);   //r_dot
    }
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
#if SUKF
  weights(0) = weights(0) + (lambda_/(lambda_+n_aug_) + 1.0-ALPHA*ALPHA+BETA);
#endif

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  S = S + R_radar_;

  //write result
  *z_out    = z_pred;
  *S_out    = S;
  *Zsig_out = Zsig;

  return true;
}

void UKF::PredictLaserMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_out)
{
    //set measurement dimension, radar can measure r, phi, and r_dot
    const int n_z = 2;

    //set vector for weights
    VectorXd weights = VectorXd(2*n_aug_+1);
#if SUKF
    const double weight_0 = lambda_/(ALPHA*ALPHA*(lambda_+n_aug_)) + 1.0/(ALPHA*ALPHA) - 1.0;
#else
    const double weight_0 = lambda_/(lambda_+n_aug_);
#endif
    weights(0) = weight_0;
    for (int i=1; i<2*n_aug_+1; i++) {  
#if SUKF
        const double weight = 0.5/(ALPHA*ALPHA*(lambda_+n_aug_));
#else
        const double weight = 0.5/(n_aug_+lambda_);
#endif
        weights(i) = weight;
    }

    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points

        // extract values for better readability
        const double p_x = Xsig_(0,i);
        const double p_y = Xsig_(1,i);

        // measurement model
        Zsig(0,i) = p_x;   //x
        Zsig(1,i) = p_y;   //y
    }

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug_+1; i++) {
        z_pred = z_pred + weights(i) * Zsig.col(i);
    }

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
#if SUKF
    weights(0) = weights(0) + (lambda_/(lambda_+n_aug_) + 1.0-ALPHA*ALPHA+BETA);
#endif
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        S = S + weights(i) * z_diff * z_diff.transpose();
    }

    //add measurement noise covariance matrix
    S = S + R_laser_;

    //write result
    *z_out    = z_pred;
    *S_out    = S;
    *Zsig_out = Zsig;
}

void UKF::InitLaser(MeasurementPackage meas_package)
{
    // state vector contains px, py, v, psi, psi_dot
    x_ = VectorXd(n_x_);
    x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0.0, 0.0, 0.0; 

    P_ << std_laspx_*std_laspx_,0,                     0,          0,            0,          
          0,                    std_laspy_*std_laspy_, 0,          0,            0,          
          0,                    0,                     V_INIT_VAR, 0,            0,          
          0,                    0,                     0,          PSI_INIT_VAR, 0,          
          0,                    0,                     0,          0,            PSI_DOT_VAR;
}
