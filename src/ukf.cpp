#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

#define EPS 0.001

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;



/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.7;

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


  // The state dimension
  n_x_ = 5;

  // The augmented state dimension
  n_aug_ = n_x_+ 2;

  // Spreading parameter
  lambda_ = 3 - n_aug_;

  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);

  // Set the predicted sigma points matrix dimentions
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Set radar measurement dimension, radar can measure r, phi, and r_dot
  n_z_R_ = 3;

  // Set lidar measurement dimension, lidar can measure x,y
  n_z_L_ = 2;

  // Add measurement noise covariance matrix
  R_Radar_ = MatrixXd(n_z_R_,n_z_R_);
  R_Radar_ <<    std_radr_*std_radr_, 0, 0,
  		         0, std_radphi_*std_radphi_, 0,
  		         0, 0,std_radrd_*std_radrd_;

  R_Lidar_ = MatrixXd(n_z_L_, n_z_L_);
  R_Lidar_ << std_laspx_*std_laspx_,0,
              0,std_laspy_*std_laspy_;
}

UKF::~UKF() {}

/**
 * @param {angle} angle The angle in rads to be normalized
 */
void UKF::AngleNorm(double *angle){
  //while ( *angle > M_PI ) *angle-=2.*M_PI;
  //while ( *angle < -M_PI ) *angle+=2.*M_PI;
  if(*angle > M_PI){
	  *angle = (int(*angle - M_PI)%int(2*M_PI)) - M_PI;
  }
  if(*angle < -M_PI){
	  *angle = (int(*angle + M_PI)%int(2*M_PI)) + M_PI;
  }
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /*****************************************************************************
  *  Initialization
  ****************************************************************************/
  if (!is_initialized_) {
	  /**
	  * Initialize the state ekf_.x_ with the first measurement.
	  * Create the covariance matrix.
	  * Convert radar from polar to cartesian coordinates.
	  */

	  // first measurement
	  cout << "EKF: " << endl;
	  // Initialize state vector
	  //x_ << 0, 0, 0, 0, 0;

	  // intialized state covariance matrix P
	  P_ << 0.5, 0, 0, 0, 0,
	        0, 0.5, 0, 0, 0,
	        0, 0, 100, 0, 0,
	        0, 0, 0, 10, 0,
	        0, 0, 0, 0, 1;

	  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		  /**
		  Convert radar from polar to cartesian coordinates and initialize state.
		  */
		  float rho = meas_package.raw_measurements_[0];
		  float phi = meas_package.raw_measurements_[1];
		  float rho_dot = meas_package.raw_measurements_[2];
		  //float px = sqrt(rho*rho/(1+tan(phi)*tan(phi)));
		  //float py = px*tan(phi);
		  float px = rho * cos(phi);
		  float py = rho * sin(phi);
		  float vx = rho_dot * cos(phi);
		  float vy = rho_dot * sin(phi);
		  float v  = sqrt(vx * vx + vy * vy);

		  x_ << px, py, v, 0, 0;

	  }
	  else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
		  /**
		  Initialize state.
		  */
		  x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
	  }

	  // Initialisation problems when px or py very close to zero
	  if (fabs(x_(0)) < EPS and fabs(x_(1)) < EPS){
		  x_(0) = EPS;
		  x_(1) = EPS;
	  }

	  // Set weights
	  double weight_0 = lambda_/(lambda_+n_aug_);
	  weights_(0) = weight_0;
	  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
	  	double weight = 0.5/(n_aug_+lambda_);
	  	weights_(i) = weight;
	  }


	  time_us_ = meas_package.timestamp_;
	  // done initializing, no need to predict or update
	  is_initialized_ = true;
	  return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/


  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
  // Save current timestamp
  time_us_ = meas_package.timestamp_;
  // Predict state and the state covariance matrix
  //Prediction(dt);
  // If timestamps are too large, break them up
  while (dt > 0.1) {
	   const double dtt = 0.1;
	   Prediction(dtt);
	   dt -= dtt;
  }
  Prediction(dt);
  //cout << "x_ = " << x_ << endl;
  //cout << "P_ = " << P_ << endl;
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
	 * Use the sensor type to perform the update step.
	 * Update the state and covariance matrices.
   */

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
	// Radar updates
	UpdateRadar(meas_package);
	//cout << "x_ = " << x_ << endl;
	//cout << "P_ = " << P_ << endl;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
	// Laser updates
	UpdateLidar(meas_package);
    // print the output
    //cout << "x_ = " << x_ << endl;
    //cout << "P_ = " << P_ << endl;
  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  // Create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  // Create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // Create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.fill(0.0);

  // Create augmented mean state
  x_aug.fill(0.0);
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // Create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  // Create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // Create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
	Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
	Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

  // Predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
	  // Extract values for better readability
	  double p_x = Xsig_aug(0,i);
	  double p_y = Xsig_aug(1,i);
	  double v = Xsig_aug(2,i);
	  double yaw = Xsig_aug(3,i);
	  double yawd = Xsig_aug(4,i);
	  double nu_a = Xsig_aug(5,i);
	  double nu_yawdd = Xsig_aug(6,i);


	  // Predicted state values
	  double px_p, py_p;

	  // Avoid division by zero
	  if (fabs(yawd) > 0.001) {
		  px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
		  py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
	  }
	  else {
		  px_p = p_x + v*delta_t*cos(yaw);
		  py_p = p_y + v*delta_t*sin(yaw);
	  }

	  double v_p = v;
	  double yaw_p = yaw + yawd*delta_t;
	  double yawd_p = yawd;

	  // Add noise
	  px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
	  py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
	  v_p = v_p + nu_a*delta_t;

	  yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
	  yawd_p = yawd_p + nu_yawdd*delta_t;

	  // Write predicted sigma point into right column
	  Xsig_pred_(0,i) = px_p;
	  Xsig_pred_(1,i) = py_p;
	  Xsig_pred_(2,i) = v_p;
	  Xsig_pred_(3,i) = yaw_p;
	  Xsig_pred_(4,i) = yawd_p;
  }

  // Predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
	  x_ = x_+ weights_(i) * Xsig_pred_.col(i);
  }

  // Predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

	  // State difference
	  VectorXd x_diff = Xsig_pred_.col(i) - x_;
	  // Angle normalization
	  AngleNorm(&(x_diff(3)));
	  P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  // Create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_L_, 2 * n_aug_ + 1);

  // Transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

      // Extract values for better readibility
	  double p_x = Xsig_pred_(0,i);
  	  double p_y = Xsig_pred_(1,i);

  	  // Measurement model
  	  Zsig(0,i) = p_x;                        //p_x
  	  Zsig(1,i) = p_y;                        //p_y
  }

  // Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_L_);

  z_pred = Zsig * weights_;

  // Measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_L_,n_z_L_);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }


  	S = S + R_Lidar_;

  	// Create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z_L_);

	// Calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

		// Residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		// State difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;

		// Angle normalization
		AngleNorm(&(x_diff(3)));

		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	 }

	 // Kalman gain K;
	 MatrixXd K = Tc * S.inverse();

	 // Measurements
	 VectorXd z = meas_package.raw_measurements_;

	 //residual
	 VectorXd z_diff = z - z_pred;

	 // Update state mean and covariance matrix
	 x_ = x_ + K * z_diff;
	 P_ = P_ - K*S*K.transpose();

	 // Calculate radar NIS
	 NIS_Lidar_ = z_diff.transpose() * S.inverse() * z_diff;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  // Create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_R_, 2 * n_aug_ + 1);

  // Transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

	  // Extract values for better readibility
	  double p_x = Xsig_pred_(0,i);
	  double p_y = Xsig_pred_(1,i);
	  double v  = Xsig_pred_(2,i);
	  double yaw = Xsig_pred_(3,i);

	  double v1 = cos(yaw)*v;
	  double v2 = sin(yaw)*v;

	  // Measurement model
	  Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
	  Zsig(1,i) = atan2(p_y,p_x);                                 //phi
	  Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
	}

	// Mean predicted measurement
	VectorXd z_pred = VectorXd(n_z_R_);
	z_pred.fill(0.0);
	for (int i=0; i < 2*n_aug_+1; i++) {
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	// Measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z_R_,n_z_R_);
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
	  //residual
	  VectorXd z_diff = Zsig.col(i) - z_pred;

	  // Angle normalization
	  AngleNorm(&(z_diff(1)));

	  S = S + weights_(i) * z_diff * z_diff.transpose();
	}


	S = S + R_Radar_;

	// Create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z_R_);

	// Calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

		// Residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		// Angle normalization
		AngleNorm(&(z_diff(1)));

	    // State difference
	    VectorXd x_diff = Xsig_pred_.col(i) - x_;

	    // Angle normalization
	    AngleNorm(&(x_diff(3)));

	    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	 }

	 // Kalman gain K;
	 MatrixXd K = Tc * S.inverse();

	 // Measurements
	 VectorXd z = meas_package.raw_measurements_;

	 //residual
	 VectorXd z_diff = z - z_pred;

	 // Angle normalization
	 AngleNorm(&(z_diff(1)));

	 // Update state mean and covariance matrix
	 x_ = x_ + K * z_diff;
	 P_ = P_ - K*S*K.transpose();

     // Calculate radar NIS
	 NIS_Radar_ = z_diff.transpose() * S.inverse() * z_diff;

}
