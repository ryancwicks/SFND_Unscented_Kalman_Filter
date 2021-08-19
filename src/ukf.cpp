#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * DO NOT MODIFY measurement noise values below.
 * These are provided by the sensor manufacturer.
 */
// Laser measurement noise standard deviation position1 in m
const double UKF::std_laspx_ = 0.15;
// Laser measurement noise standard deviation position2 in m
const double UKF::std_laspy_ = 0.15;
// Radar measurement noise standard deviation radius in m
const double UKF::std_radr_ = 0.3;
// Radar measurement noise standard deviation angle in rad
const double UKF::std_radphi_ = 0.03;
// Radar measurement noise standard deviation radius change in m/s
const double UKF::std_radrd_ = 0.3;

/**
 * End DO NOT MODIFY section for measurement noise values 
 */
  

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() :
  is_initialized_(false),
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ (true),
  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ (true),
  n_x_(5),
  n_aug_(7),
  // initial state vector
  x_ (VectorXd(n_x_)),
  // initial covariance matrix
  P_ (MatrixXd(n_x_, n_x_)),
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ (30),
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ (30),
  //weights
  weights_(VectorXd(2*n_aug_+1)),
  lambda_(3 - n_aug_)
{
  // set weights
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < 2*n_aug_+1; ++i){
    weights_[i] = 1./2./(lambda_+n_aug_);
  }
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(const MeasurementPackage& meas_package) {
  if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER ) {
    UpdateLidar(meas_package);
  } else if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }
}

void UKF::Prediction(double delta_t) {
  if (is_initialized_) {
    MatrixXd Xsig_aug;
    generateAugmentedSigmaPointMatrix(Xsig_aug);
    generateSigmaPointPrediction(Xsig_aug, delta_t);
    generatePredictedMeanAndCovariance();
  }
}


void UKF::generateSigmaPointMatrix( MatrixXd & Xsig) {
  // create sigma point matrix
  Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);
  lambda_ = 3 - n_x_;
  
  // calculate square root of P
  MatrixXd A = P_.llt().matrixL();

  // calculate sigma points ...
  MatrixXd x_5d = MatrixXd(n_x_, n_x_);
  for (int i = 0; i < n_x_; ++i) {
      x_5d.col(i) = x_;
  }
  Xsig.col(0) = x_;
  Xsig.block(0, 1, n_x_, n_x_) = x_5d + std::sqrt(lambda_+n_x_) * A;
  Xsig.block(0, n_x_+1, n_x_, n_x_) = x_5d - std::sqrt(lambda_+n_x_)*A;
}

void UKF::generateAugmentedSigmaPointMatrix( MatrixXd & Xsig_aug) {
  // create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  
  //Update lambda
  lambda_ = 3 - n_aug_;

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // create sigma point matrix
  Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
 
  // create augmented mean state
  x_aug.topLeftCorner(n_x_, 1) =  x_;
  x_aug(n_x_) = 0.0;
  x_aug(n_x_+1) = 0.0;
  
  MatrixXd x_aug_block = MatrixXd(n_aug_, n_aug_);
  for (int i = 0; i < n_aug_; ++i) {
      x_aug_block.col(i) = x_aug;
  }

  // create augmented covariance matrix
  P_aug = MatrixXd::Zero(n_aug_, n_aug_);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_aug_-2, n_aug_-2) = std_a_*std_a_;
  P_aug(n_aug_-1, n_aug_-1) = std_yawdd_*std_yawdd_;

  // create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  // create augmented sigma points
  Xsig_aug.col(0) =  x_aug;
  Xsig_aug.block(0, 1, n_aug_, n_aug_) = x_aug_block + std::sqrt(lambda_ + n_aug_)*A;
  Xsig_aug.block(0, n_aug_+1, n_aug_, n_aug_) = x_aug_block - std::sqrt(lambda_ + n_aug_)*A;

}

void UKF::generateSigmaPointPrediction ( const Eigen::MatrixXd & Xsig_aug, double delta_t) {

  // create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // predict sigma points
  for (int i = 0; i < 2*n_aug_+1; ++i) {
      double px = Xsig_aug.col(i)(0);
      double py = Xsig_aug.col(i)(1);
      double v = Xsig_aug.col(i)(2);
      double psi = Xsig_aug.col(i)(3);
      double psi_d = Xsig_aug.col(i)(4);
      double mu_a = Xsig_aug.col(i)(5);
      double mu_psi_dd = Xsig_aug.col(i)(6);
      
      Xsig_pred_.col(i) << px, py, v, psi, psi_d;
      
      if (std::abs(psi_d) > 0.00000001) { //prevent division by 0.
          Xsig_pred_.col(i)(0) += v/psi_d * (sin(psi+psi_d*delta_t)-sin(psi))  + 1.0/2.0*delta_t*delta_t * cos(psi)*mu_a;
          Xsig_pred_.col(i)(1) += v/psi_d * (-cos(psi+psi_d*delta_t)+cos(psi)) + 1.0/2.0*delta_t*delta_t * sin(psi)*mu_a;;
          Xsig_pred_.col(i)(2) += delta_t * mu_a;
          Xsig_pred_.col(i)(3) += psi_d * delta_t + 1.0/2.0*delta_t*delta_t*mu_psi_dd;
          Xsig_pred_.col(i)(4) += delta_t * mu_psi_dd;
      } else {
          Xsig_pred_.col(i)(0) += v * cos(psi)*delta_t  + 1.0/2.0*delta_t*delta_t * cos(psi)*mu_a;
          Xsig_pred_.col(i)(1) += v * sin(psi)*delta_t  + 1.0/2.0*delta_t*delta_t * sin(psi)*mu_a;;
          Xsig_pred_.col(i)(2) += delta_t * mu_a;
          Xsig_pred_.col(i)(3) += psi_d * delta_t + 1.0/2.0*delta_t*delta_t*mu_psi_dd;
          Xsig_pred_.col(i)(4) += delta_t * mu_psi_dd;
      }
  }

}

void UKF::generatePredictedMeanAndCovariance() {
  x_ = VectorXd::Zero(n_x_);
  P_ = MatrixXd::Zero(n_x_, n_x_);

    for (int i = 0; i < 2*n_aug_+1; ++i) {
        x_ += weights_[i] * Xsig_pred_.col(i);
    }
    for (int i = 0; i < 2*n_aug_+1; ++i) {
        P_ += weights_[i] * (Xsig_pred_.col(i) - x_) * (Xsig_pred_.col(i) - x_).transpose();
    }
}

void UKF::UpdateLidar(const MeasurementPackage& meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
}

void UKF::predictLidarMeasurement( Eigen::VectorXd & z_out, Eigen::MatrixXd & S_out) {

}

void UKF::UpdateRadar(const MeasurementPackage& meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}

void UKF::predictRadarMeasurement( Eigen::VectorXd & z_out, Eigen::MatrixXd & S_out) {

}