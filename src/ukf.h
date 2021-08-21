#ifndef UKF_H
#define UKF_H

#include "Eigen/Dense"
#include "measurement_package.h"

class UKF {
 public:
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
   */
  void ProcessMeasurement(const MeasurementPackage& meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);


  // initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  // if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  // if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  // State dimension
  int n_x_;

  // Augmented state dimension
  int n_aug_;

  // measurement dimension radar
  int n_z_radar_;

  // measurement dimension lidar
  int n_z_lidar_;

  // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // predicted sigma points matrix
  Eigen::MatrixXd Xsig_pred_;

  // time when the state is true, in us
  long long time_us_;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  // Laser measurement noise standard deviation position1 in m
  const static double std_laspx_;

  // Laser measurement noise standard deviation position2 in m
  const static double std_laspy_;

  // Radar measurement noise standard deviation radius in m
  const static double std_radr_;

  // Radar measurement noise standard deviation angle in rad
  const static double std_radphi_;

  // Radar measurement noise standard deviation radius change in m/s
  const static double std_radrd_ ;

  // Weights of sigma points
  Eigen::VectorXd weights_;

  // Sigma point spreading parameter
  double lambda_;

  //Measurement matrix for the laser.
  Eigen::MatrixXd H_laser_; 

  protected:

    /**
     * Function to generate sigma points (not actually used, kept for completeness).
     * @param Xsig [in/out] sigma point matrix, returned by reference
     */
    void generateSigmaPointMatrix( Eigen::MatrixXd & Xsig);

    /**
     * Function to generate augmented sigma points .
     * @param Xsig_aug [in/out] augmented sigma point matrix, returned by reference
     */
    void generateAugmentedSigmaPointMatrix( Eigen::MatrixXd & Xsig_aug);

    /**
     * Predict sigma points for update step. Stores the result in the class Xsig_pred_ variable.
     * @param Xsig_aug [in] augmented sigma point matrix
     * @param delta_t [in] time step in seconds
     */
    void generateSigmaPointPrediction ( const Eigen::MatrixXd & Xsig_aug, double delta_t);

    /**
     * Generate the new mean and covariance from the sigma point prediction. This function relies entirely on the internal state.
     */
    void generatePredictedMeanAndCovariance();

    /**
     * Updates the state and the state covariance matrix using a laser measurement with 
     * a UKF update. Not used in this project.
     * @param meas_package The measurement at k+1
     */
    void UpdateLidarUKF(const MeasurementPackage& meas_package);

    /**
     * Updates the state and the state covariance matrix using a laser measurement with a 
     * linear filter update.
     * @param meas_package The measurement at k+1
     */
    void UpdateLidarLinear(const MeasurementPackage& meas_package);

    /**
     * Updates the state and the state covariance matrix using a radar measurement
     * @param meas_package The measurement at k+1
     */
    void UpdateRadar(const MeasurementPackage& meas_package);

    /**
     * Get the predicted radar measurement state.
     * @param z_out [out] predicted radar measurement mean
     * @param S_out [out] predicted radar covariance state
     * @param Z_sig [out] sigma points in measurement space
     */
    void predictRadarMeasurement( Eigen::VectorXd & z_out, Eigen::MatrixXd & S_out, Eigen::MatrixXd & Z_sig);

    /**
     * Get the predicted LIDAR measurement state. UKF model. Not used in this project.
     * @param z_out [out] predicted lidar measurement mean
     * @param X_out [out] predicted lidar covariance state
     */
    void predictLidarMeasurementUKF( Eigen::VectorXd & z_out, Eigen::MatrixXd & S_out, Eigen::MatrixXd & Z_sig);

    /**
     * Get the predicted LIDAR measurement state. Linear model.
     * @param z_out [out] predicted lidar measurement mean
     * @param X_out [out] predicted lidar covariance state
     */
    void predictLidarMeasurementLinear( Eigen::VectorXd & z_out, Eigen::MatrixXd & S_out);


};

#endif  // UKF_H