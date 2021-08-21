#include "UKFTest.h"
#include "measurement_package.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

const double UKFTest::SMALL_VALUE = 0.00001;

void UKFTest::testGenerateSigmaPointMatrix () {

    std::cout << "Testing Sigma Point Generation..." << std::endl;
    // set example state
    x_ = VectorXd(n_x_);
    x_ <<   5.7441,
            1.3800,
            2.2049,
            0.5015,
            0.3528;

    // set example covariance matrix
    P_ = MatrixXd(n_x_, n_x_);
    P_ <<    0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
            -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
            0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
            -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
            -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

    // Expected Result
    MatrixXd Xsig_expected = MatrixXd(n_x_, 2 * n_x_ + 1);
    Xsig_expected << 5.7441,  5.85768,   5.7441,   5.7441,  5.7441,   5.7441,  5.63052,   5.7441,   5.7441,   5.7441,   5.7441,
                     1.38,  1.34566,  1.52806,    1.38,     1.38,     1.38,   1.41434,  1.23194,     1.38 ,    1.38,     1.38,
                     2.2049,  2.28414,  2.24557,  2.29582,   2.2049,   2.2049,  2.12566,  2.16423,  2.11398,   2.2049,   2.2049,
                     0.5015,  0.44339, 0.631886, 0.516923, 0.595227,   0.5015,  0.55961, 0.371114, 0.486077, 0.407773,   0.5015,
                     0.3528, 0.299973, 0.462123, 0.376339,  0.48417, 0.418721, 0.405627, 0.243477, 0.329261,  0.22143, 0.286879;

    MatrixXd Xsig;

    generateSigmaPointMatrix(Xsig);

    assert (Xsig.size() == Xsig_expected.size());
    assert (Xsig.cols() == Xsig_expected.cols());
    assert (Xsig.rows() == Xsig_expected.rows());
    for (auto row = 0; row < Xsig.rows(); ++row) {
        for (auto col = 0; col < Xsig.cols(); ++col) {
            assert(almostEqual(Xsig(row, col), Xsig_expected(row, col)));
        }
    }

    std::cout << "Sigma Point generation test passed." << std::endl;
}

void UKFTest::testGenerateAugmentedSigmaPointMatrix() {
    
    std::cout << "Testing Augmented Sigma Point Generation..." << std::endl;

    // set example state
    VectorXd x_ = VectorXd(n_x_);
    x_ <<   5.7441,
            1.3800,
            2.2049,
            0.5015,
            0.3528;

    // create example covariance matrix
    P_ = MatrixXd(n_x_, n_x_);
    P_ <<   0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
            -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
            0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
            -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
            -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

    std_a_ = 0.2;
    std_yawdd_ = 0.2;

    MatrixXd Xsig_aug_expt = MatrixXd(n_aug_, 2*n_aug_ + 1);
    Xsig_aug_expt << 5.7441,  5.85768,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,  5.63052,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,
                       1.38,  1.34566,  1.52806,     1.38,     1.38,     1.38,     1.38,     1.38,  1.41434,  1.23194,     1.38,     1.38,     1.38,     1.38,     1.38,
                     2.2049,  2.28414,  2.24557,  2.29582,   2.2049,   2.2049,   2.2049,   2.2049,  2.12566,  2.16423,  2.11398,   2.2049,   2.2049,   2.2049,   2.2049,
                     0.5015,  0.44339, 0.631886, 0.516923, 0.595227,   0.5015,   0.5015,   0.5015,  0.55961, 0.371114, 0.486077, 0.407773,   0.5015,   0.5015,   0.5015,
                     0.3528, 0.299973, 0.462123, 0.376339,  0.48417, 0.418721,   0.3528,   0.3528, 0.405627, 0.243477, 0.329261,  0.22143, 0.286879,   0.3528,   0.3528,
                          0,        0,        0,        0,        0,        0,  0.34641,        0,        0,        0,        0,        0,        0, -0.34641,        0,
                          0,        0,        0,        0,        0,        0,        0,  0.34641,        0,        0,        0,        0,        0,        0, -0.34641;
 
    MatrixXd Xsig_aug;
    generateAugmentedSigmaPointMatrix(Xsig_aug);

    assert (Xsig_aug.size() == Xsig_aug_expt.size());
    assert (Xsig_aug.cols() == Xsig_aug_expt.cols());
    assert (Xsig_aug.rows() == Xsig_aug_expt.rows());
    for (auto row = 0; row < Xsig_aug.rows(); ++row) {
        for (auto col = 0; col < Xsig_aug.cols(); ++col) {
            assert(almostEqual(Xsig_aug(row, col), Xsig_aug_expt(row, col)));
        }
    }

    std::cout << "Augmented Sigma Point generation test passed." << std::endl;
}

void UKFTest::testGenerateSigmaPointPrediction() {
    std::cout << "Testing Augmented Sigma Point Prediction..." << std::endl;

    // set state dimension
    n_x_ = 5;

    // set augmented dimension
    n_aug_ = 7;

    // set lambda
    lambda_ = 3 - n_aug_;

    // create example sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    Xsig_aug <<
      5.7441,  5.85768,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.63052,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,
        1.38,  1.34566,  1.52806,     1.38,     1.38,     1.38,     1.38,     1.38,   1.41434,  1.23194,     1.38,     1.38,     1.38,     1.38,     1.38,
      2.2049,  2.28414,  2.24557,  2.29582,   2.2049,   2.2049,   2.2049,   2.2049,   2.12566,  2.16423,  2.11398,   2.2049,   2.2049,   2.2049,   2.2049,
      0.5015,  0.44339, 0.631886, 0.516923, 0.595227,   0.5015,   0.5015,   0.5015,   0.55961, 0.371114, 0.486077, 0.407773,   0.5015,   0.5015,   0.5015,
      0.3528, 0.299973, 0.462123, 0.376339,  0.48417, 0.418721,   0.3528,   0.3528,  0.405627, 0.243477, 0.329261,  0.22143, 0.286879,   0.3528,   0.3528,
           0,        0,        0,        0,        0,        0,  0.34641,        0,         0,        0,        0,        0,        0, -0.34641,        0,
           0,        0,        0,        0,        0,        0,        0,  0.34641,         0,        0,        0,        0,        0,        0, -0.34641;


    MatrixXd Xsig_pred_expt = MatrixXd(n_x_, 2 * n_aug_ + 1);
    Xsig_pred_expt << 5.93553, 6.06251,   5.92217, 5.9415, 5.92361, 5.93516, 5.93705, 5.93553, 5.80832, 5.94481, 5.92935, 5.94553, 5.93589, 5.93401, 5.93553,
                      1.48939, 1.44673,   1.66484, 1.49719, 1.508, 1.49001, 1.49022, 1.48939, 1.5308, 1.31287, 1.48182, 1.46967, 1.48876, 1.48855, 1.48939,
                       2.2049, 2.28414,   2.24557, 2.29582,  2.2049, 2.2049, 2.23954, 2.2049, 2.12566, 2.16423, 2.11398, 2.2049, 2.2049, 2.17026, 2.2049,
                      0.53678, 0.473387, 0.678098, 0.554557, 0.643644, 0.543372, 0.53678, 0.538512, 0.600173, 0.395462, 0.519003, 0.429916, 0.530188, 0.53678, 0.535048,
                       0.3528, 0.299973, 0.462123, 0.376339, 0.48417, 0.418721, 0.3528, 0.387441, 0.405627, 0.243477, 0.329261, 0.22143, 0.286879, 0.3528, 0.318159;


    double delta_t = 0.1; // time diff in sec
    generateSigmaPointPrediction(Xsig_aug, delta_t);

    assert (Xsig_pred_.size() == Xsig_pred_expt.size());
    assert (Xsig_pred_.cols() == Xsig_pred_expt.cols());
    assert (Xsig_pred_.rows() == Xsig_pred_expt.rows());
    for (auto row = 0; row < Xsig_pred_.rows(); ++row) {
        for (auto col = 0; col < Xsig_pred_.cols(); ++col) {
            assert(almostEqual(Xsig_pred_(row, col), Xsig_pred_expt(row, col)));
        }
    }

    std::cout << "Augmented Sigma Point generation test passed." << std::endl;
}

void UKFTest::testPredictMeanAndCovariance() {
    std::cout << "Testing mean and covariance prediction..." << std::endl;

    // create example matrix with predicted sigma points
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    Xsig_pred_ <<
         5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
           1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
          2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
         0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
          0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

    VectorXd x_expt = VectorXd(n_x_);
    x_expt <<  5.93637, 1.49035, 2.20528, 0.536853, 0.35357;
    MatrixXd P_expt = MatrixXd(n_x_, n_x_);
    P_expt <<  0.00543425, -0.0024053, 0.00341576, -0.00348196, -0.00299378,
               -0.0024053, 0.010845, 0.0014923, 0.00980182, 0.00791091,
               0.00341576, 0.0014923, 0.00580129, 0.000778632, 0.000792973,
               -0.00348196, 0.00980182, 0.000778632, 0.0119238, 0.0112491,
               -0.00299378, 0.00791091, 0.000792973, 0.0112491, 0.0126972;


    generatePredictedMeanAndCovariance();

    assert (x_.size() == x_expt.size());
    for (auto i = 0; i < x_.rows(); ++i) {
        assert (almostEqual(x_(i), x_expt(i)));
    }
    assert (P_.size() == P_expt.size());
    assert (P_.cols() == P_expt.cols());
    assert (P_.rows() == P_expt.rows());
    for (auto row = 0; row < P_.rows(); ++row) {
        for (auto col = 0; col < P_.cols(); ++col) {
            assert(almostEqual(P_(row, col), P_expt(row, col)));
        }
    }

    std::cout << "Predict mean and covariance test passed." << std::endl;
}

void UKFTest::testRadarMeasurementPrediction() {
    std::cout << "Testing Radar measurement prediction..." << std::endl;

    // create example matrix with predicted sigma points
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    Xsig_pred_ <<
         5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
           1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
          2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
         0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
          0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

    VectorXd z_expt = VectorXd(n_z_radar_);
    z_expt <<  6.12155, 0.245993, 2.10313;
    MatrixXd S_expt (n_z_radar_, n_z_radar_);
    S_expt <<   0.0946171, -0.000139448, 0.00407016,
                -0.000139448, 0.0012113, -0.000770652,
                0.00407016, -0.000770652, 0.0980917;

    VectorXd z_out;
    MatrixXd S_out;
    MatrixXd Z_sig;

    predictRadarMeasurement(z_out, S_out, Z_sig);

    assert (z_out.size() == z_expt.size());
    for (auto i = 0; i < z_out.rows(); ++i) {
        assert (almostEqual(z_out(i), z_expt(i)));
    }
    assert (S_out.size() == S_expt.size());
    assert (S_out.cols() == S_expt.cols());
    assert (S_out.rows() == S_expt.rows());
    for (auto row = 0; row < S_out.rows(); ++row) {
        for (auto col = 0; col < S_out.cols(); ++col) {
            assert(almostEqual(S_out(row, col), S_expt(row, col)));
        }
    }

    std::cout << "Radar measurement prediction test passed" << std::endl;
}

void UKFTest::testRadarMeasurementUpdate() {

    std::cout << "Testing Radar measurement prediction..." << std::endl;

    // create example matrix with predicted sigma points in state space
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    Xsig_pred_ <<
     5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
       1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
      2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
     0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
      0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

    // create example vector for predicted state mean
    x_ = VectorXd(n_x_);
    x_ <<
     5.93637,
     1.49035,
     2.20528,
    0.536853,
    0.353577;

    // create example matrix for predicted state covariance
    P_ = MatrixXd(n_x_,n_x_);
    P_ <<
    0.0054342,  -0.002405,  0.0034157, -0.0034819, -0.00299378,
    -0.002405,    0.01084,   0.001492,  0.0098018,  0.00791091,
    0.0034157,   0.001492,  0.0058012, 0.00077863, 0.000792973,
   -0.0034819,  0.0098018, 0.00077863,   0.011923,   0.0112491,
   -0.0029937,  0.0079109, 0.00079297,   0.011249,   0.0126972;

    VectorXd z = VectorXd(n_z_radar_);
    z << 5.9214,   // rho in m
         0.2187,   // phi in rad
         2.0062;   // rho_dot in m/s

    MeasurementPackage measure {0, MeasurementPackage::RADAR,z};

    /**VectorXd x_expt = VectorXd(n_x_);
    x_expt <<  5.92276, 1.41823, 2.15593, 0.489274, 0.321338;
    MatrixXd P_expt = MatrixXd(n_x_, n_x_);
    P_expt <<  0.00361579, -0.000357881, 0.00208316, -0.000937196, -0.00071727,
               -0.000357881, 0.00539867, 0.00156846, 0.00455342, 0.00358885,
               0.00208316, 0.00156846, 0.00410651, 0.00160333, 0.00171811,
               -0.000937196, 0.00455342, 0.00160333, 0.00652634, 0.00669436,
               -0.00071719, 0.00358884, 0.00171811, 0.00669426, 0.00881797;
    */
    UpdateRadar(measure);

    /**
    std::cout << x_ << std::endl;
    std::cout << P_ << std::endl;


    assert (x_.size() == x_expt.size());
    for (auto i = 0; i < x_.rows(); ++i) {
        assert (almostEqual(x_(i), x_expt(i)));
    }
    assert (P_.size() == P_expt.size());
    assert (P_.cols() == P_expt.cols());
    assert (P_.rows() == P_expt.rows());
    for (auto row = 0; row < P_.rows(); ++row) {
        for (auto col = 0; col < P_.cols(); ++col) {
            assert(almostEqual(P_(row, col), P_expt(row, col)));
        }
    }
    */

    std::cout << "Radar measurement update test NOT IMPLEMENTED" << std::endl;
}

void UKFTest::testLidarMeasurementPrediction() {
    std::cout << "Testing Lidar measurement prediction..." << std::endl;

    x_ = VectorXd(n_x_);
    x_ <<  5.93637, 1.49035, 2.20528, 0.536853, 0.35357;
    P_ = MatrixXd(n_x_, n_x_);
    P_ <<  0.00543425, -0.0024053, 0.00341576, -0.00348196, -0.00299378,
               -0.0024053, 0.010845, 0.0014923, 0.00980182, 0.00791091,
               0.00341576, 0.0014923, 0.00580129, 0.000778632, 0.000792973,
               -0.00348196, 0.00980182, 0.000778632, 0.0119238, 0.0112491,
               -0.00299378, 0.00791091, 0.000792973, 0.0112491, 0.0126972;

    VectorXd z_expt = VectorXd(n_z_lidar_);
    z_expt <<  5.93637, 1.49035;
    MatrixXd S_expt (n_z_lidar_, n_z_lidar_);
    S_expt << 0.00543425+std_laspx_*std_laspx_, -0.0024053,
              -0.0024053, 0.010845 + std_laspy_*std_laspy_; 
              

    VectorXd z_out;
    MatrixXd S_out;
    MatrixXd Z_sig;

    predictLidarMeasurementLinear(z_out, S_out);
    //predictLidarMeasurementUKF(z_out, S_out, Z_sig);

    assert (z_out.size() == z_expt.size());
    for (auto i = 0; i < z_out.rows(); ++i) {
        assert (almostEqual(z_out(i), z_expt(i)));
    }
    assert (S_out.size() == S_expt.size());
    assert (S_out.cols() == S_expt.cols());
    assert (S_out.rows() == S_expt.rows());
    for (auto row = 0; row < S_out.rows(); ++row) {
        for (auto col = 0; col < S_out.cols(); ++col) {
            assert(almostEqual(S_out(row, col), S_expt(row, col)));
        }
    }

    std::cout << "Lidar measurement prediction test passed" << std::endl;
}

void UKFTest::testLidarMeasurementUpdate() {
    std::cout << "Testing Lidar measurement update..." << std::endl;

    // create example vector for predicted state mean
    x_ = VectorXd(n_x_);
    x_ <<
     5.93637,
     1.49035,
     2.20528,
    0.536853,
    0.353577;

    // create example matrix for predicted state covariance
    P_ = MatrixXd(n_x_,n_x_);
    P_ <<
    0.0054342,  -0.002405,  0.0034157, -0.0034819, -0.00299378,
    -0.002405,    0.01084,   0.001492,  0.0098018,  0.00791091,
    0.0034157,   0.001492,  0.0058012, 0.00077863, 0.000792973,
   -0.0034819,  0.0098018, 0.00077863,   0.011923,   0.0112491,
   -0.0029937,  0.0079109, 0.00079297,   0.011249,   0.0126972;

    VectorXd z = VectorXd(n_z_lidar_);
    z << 6.5, 1.9;

    MeasurementPackage measure {0, MeasurementPackage::LASER, z};

    /**VectorXd x_expt = VectorXd(n_x_);
    x_expt <<  5.92276, 1.41823, 2.15593, 0.489274, 0.321338;
    MatrixXd P_expt = MatrixXd(n_x_, n_x_);
    P_expt <<  0.00361579, -0.000357881, 0.00208316, -0.000937196, -0.00071727,
               -0.000357881, 0.00539867, 0.00156846, 0.00455342, 0.00358885,
               0.00208316, 0.00156846, 0.00410651, 0.00160333, 0.00171811,
               -0.000937196, 0.00455342, 0.00160333, 0.00652634, 0.00669436,
               -0.00071719, 0.00358884, 0.00171811, 0.00669426, 0.00881797;
    */
    UpdateLidarLinear(measure);

    std::cout << "Lidar measurement update test NOT IMPLEMENTED" << std::endl;
}


bool UKFTest::almostEqual (double a, double b) const {
    return std::abs(a - b) < SMALL_VALUE;
}