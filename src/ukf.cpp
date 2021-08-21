#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

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
  n_z_radar_(3),
  n_z_lidar_(2),
  // initial state vector
  x_ (VectorXd(n_x_)),
  // initial covariance matrix
  P_ (MatrixXd(n_x_, n_x_)),
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ (2.0),
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ (1.0),
  //weights
  weights_(VectorXd(2*n_aug_+1)),
  lambda_(3.0 - n_aug_)
{
  // set weights
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < 2*n_aug_+1; ++i){
    weights_[i] = 1./2./(lambda_+n_aug_);
  }
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(const MeasurementPackage& meas_package) {
  if (!is_initialized_) {
    time_us_ = meas_package.timestamp_;
    x_ = VectorXd::Zero(n_x_);
    if (meas_package.sensor_type_ == meas_package.LASER) {
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
    } else {
      x_(0) = meas_package.raw_measurements_(0)*std::cos(meas_package.raw_measurements_(1));
      x_(1) = meas_package.raw_measurements_(0)*std::sin(meas_package.raw_measurements_(1));
    }
    P_ = 100.0 * MatrixXd::Identity(n_x_, n_x_);
    P_(2, 2) = 25.0;
    P_(3, 3) = 1.0;
    P_(4, 4) = 1.0;
    is_initialized_ = true;  
    MatrixXd Xsig_aug;
    generateAugmentedSigmaPointMatrix(Xsig_aug);
    generateSigmaPointPrediction(Xsig_aug, 0);
    generatePredictedMeanAndCovariance();
  } 
  double diff_t = (meas_package.timestamp_ - time_us_)/1.0E6;
  Prediction(diff_t);
  time_us_ = meas_package.timestamp_;
  if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER ) {
    UpdateLidarLinear(meas_package);
  } else if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }
}

void UKF::Prediction(double delta_t) {
  MatrixXd Xsig_aug;
  generateAugmentedSigmaPointMatrix(Xsig_aug);
  generateSigmaPointPrediction(Xsig_aug, delta_t);
  generatePredictedMeanAndCovariance();

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

  lambda_ = 3.0-n_aug_;

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
  x_ = VectorXd::Zero(n_x_);
      
      Xsig_pred_.col(i) << px, py, v, psi, psi_d;
      
      if (std::abs(psi_d) > 0.001) { //prevent division by 0.
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

      //Xsig_pred_.col(i)(3) = Xsig_pred_.col(i)(3) - int64_t(Xsig_pred_.col(i)(3) / (2*M_PI)) * 2.0*M_PI; 
      //Xsig_pred_.col(i)(4) = Xsig_pred_.col(i)(4) - int64_t(Xsig_pred_.col(i)(4) / (2*M_PI)) * 2.0*M_PI;
  }

}

void UKF::generatePredictedMeanAndCovariance() {
  x_ = VectorXd::Zero(n_x_);
  P_ = MatrixXd::Zero(n_x_, n_x_);

    for (int i = 0; i < 2*n_aug_+1; ++i) {
        x_ += weights_[i] * Xsig_pred_.col(i);
    }
    for (int i = 0; i < 2*n_aug_+1; ++i) {
      // state difference
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      // angle normalization
      double yaw = x_diff(3);
      while (yaw> M_PI) yaw-=2.*M_PI;
      while (yaw<-M_PI) yaw+=2.*M_PI;
      x_diff(3) = yaw;

      P_ += weights_(i) * x_diff * x_diff.transpose() ;
    }
}

void UKF::UpdateLidarLinear(const MeasurementPackage& meas_package) {

  VectorXd z_pred;
  MatrixXd S;

  //set up the measurement matrix
  MatrixXd H = MatrixXd (n_z_lidar_, n_x_);
  H << 1.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 1.0, 0.0, 0.0, 0.0;

  predictLidarMeasurementLinear(z_pred, S);

  // calculate Kalman gain K;
  MatrixXd K;
  K = P_*H.transpose()*S.inverse();

  // update state mean and covariance matrix
  VectorXd z_diff = meas_package.raw_measurements_-z_pred;
  x_ += K*(z_diff);
  P_ -= K*H*P_;
}

void UKF::predictLidarMeasurementLinear( Eigen::VectorXd & z_out, Eigen::MatrixXd & S_out) {
  //Assumes the state has already been updated to the current measurement time.
  // mean predicted measurement
  z_out = VectorXd(n_z_lidar_);
  
  // measurement covariance matrix S
  S_out = MatrixXd(n_z_lidar_, n_z_lidar_);

  //set up the measurement matrix
  MatrixXd H = MatrixXd (n_z_lidar_, n_x_);
  H << 1.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 1.0, 0.0, 0.0, 0.0;

  //measured noise covariance matrix
  MatrixXd R = MatrixXd(n_z_lidar_, n_z_lidar_);
  R << std_laspx_*std_laspx_, 0, 
       0, std_laspy_*std_laspy_;

  z_out = H*x_;
  S_out = H*P_*H.transpose() + R;
}


void UKF::UpdateLidarUKF(const MeasurementPackage& meas_package) {

  VectorXd z_pred;
  MatrixXd S;
  MatrixXd Z_sig;

  predictLidarMeasurementUKF(z_pred, S, Z_sig);

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z_lidar_);

  // calculate cross correlation matrix
  for (int i = 0; i<2*n_aug_+1; ++i) {
    // residual
    VectorXd z_diff = Z_sig.col(i) - z_pred;
  
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc += weights_[i] * (x_diff) * (z_diff).transpose();
  }

  // calculate Kalman gain K;
  MatrixXd K;
  K = Tc*S.inverse();

  // residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

  // update state mean and covariance matrix
  x_ += K*(z_diff);
  
  P_ -= K*S*K.transpose();
}

void UKF::predictLidarMeasurementUKF( Eigen::VectorXd & z_out, Eigen::MatrixXd & S_out, Eigen::MatrixXd & Z_sig) {
  //Assumes the state has already been updated to the current measurement time.
  // mean predicted measurement
  z_out = VectorXd(n_z_lidar_);
  
  // measurement covariance matrix S
  S_out = MatrixXd(n_z_lidar_, n_z_lidar_);

  // sigma point matrix 
  Z_sig = MatrixXd(n_z_lidar_, 2 * n_aug_ + 1);

     // transform sigma points into measurement space
    for (int i = 0; i < 2*n_aug_+1; ++i) {
        double px = Xsig_pred_.col(i)[0];
        double py = Xsig_pred_.col(i)[1];
        
        Z_sig.col(i) << px, py;
    }
  
    // calculate mean predicted measurement
    z_out = VectorXd::Zero(n_z_lidar_);
    for (int i = 0; i < 2*n_aug_+1; ++i){
        z_out += weights_[i] * Z_sig.col(i);
    }
    S_out = MatrixXd::Zero(n_z_lidar_, n_z_lidar_);
    for (int i = 0; i < 2*n_aug_+1; ++i) {
        VectorXd z_diff = Z_sig.col(i) - z_out;
        
        S_out += weights_[i] * (z_diff)*(z_diff).transpose();
    }
  

  //measured noise covariance matrix
  MatrixXd R = MatrixXd(n_z_lidar_, n_z_lidar_);
  R << std_laspx_*std_laspx_, 0, 
       0, std_laspy_*std_laspy_;

  S_out += R;
}

void UKF::UpdateRadar(const MeasurementPackage& meas_package) {
  VectorXd z_pred;
  MatrixXd S;
  MatrixXd Z_sig;

  predictRadarMeasurement(z_pred, S, Z_sig);

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z_radar_);

  // calculate cross correlation matrix
  for (int i = 0; i<2*n_aug_+1; ++i) {
    // residual
    VectorXd z_diff = Z_sig.col(i) - z_pred;
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
 
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc += weights_[i] * (x_diff) * (z_diff).transpose();
  }

  // calculate Kalman gain K;
  MatrixXd K;
  K = Tc*S.inverse();

  // residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

  // angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // update state mean and covariance matrix
  x_ += K*(z_diff);
  
  P_ -= K*S*K.transpose();
}

void UKF::predictRadarMeasurement( Eigen::VectorXd & z_out, Eigen::MatrixXd & S_out, Eigen::MatrixXd & Z_sig) {
    //Assumes the state has already been updated to the current time.

    // create matrix for sigma points in measurement space
    Z_sig = MatrixXd(n_z_radar_, 2 * n_aug_ + 1);

    // mean predicted measurement
    z_out = VectorXd(n_z_radar_);
  
    // measurement covariance matrix S
    S_out = MatrixXd(n_z_radar_,n_z_radar_);

    // transform sigma points into measurement space
    for (int i = 0; i < 2*n_aug_+1; ++i) {
        double px = Xsig_pred_.col(i)[0];
        double py = Xsig_pred_.col(i)[1];
        double v = Xsig_pred_.col(i)[2];
        double psi = Xsig_pred_.col(i)[3];
        //double psi_d = Xsig_pred_.col(i)[4];
            
        double r = sqrt(px*px + py*py);
        double theta = atan2(py, px);
        double r_d(0.0);
        if (std::abs(r) > 0.00001) {
          r_d = v / r * (px*cos(psi) + py*sin(psi));
        } else {
          r_d = v;
        }
        Z_sig.col(i) << r, theta, r_d;
    }
  
    // calculate mean predicted measurement
    z_out = VectorXd::Zero(n_z_radar_);
    for (int i = 0; i < 2*n_aug_+1; ++i){
        z_out += weights_[i] * Z_sig.col(i);
    }
    S_out = MatrixXd::Zero(n_z_radar_, n_z_radar_);
    for (int i = 0; i < 2*n_aug_+1; ++i) {
        VectorXd z_diff = Z_sig.col(i) - z_out;
        
        // angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        S_out += weights_[i] * (z_diff)*(z_diff).transpose();
    }
  
    MatrixXd R = MatrixXd(n_z_radar_, n_z_radar_);
    R << std_radr_*std_radr_, 0, 0, 
         0, std_radphi_*std_radphi_, 0, 
         0, 0, std_radrd_*std_radrd_;
    S_out += R;
}