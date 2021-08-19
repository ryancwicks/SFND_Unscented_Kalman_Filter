#include "UKFTest.h"
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

bool UKFTest::almostEqual (double a, double b) const {
    return std::abs(a - b) < SMALL_VALUE;
}