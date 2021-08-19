#pragma once
/**
 * UKFTest.h
 * 
 * This class is used to test the UKF class components.
 * 
 */


#include "ukf.h"
#include "Eigen/Dense"

class UKFTest : public UKF {
public:

    /**
     * Run the test for the Sigma Point generation.
     */
    void testGenerateSigmaPointMatrix ();

    /**
     * Run the test for the Augmented Sigma Point generation.
     */
    void testGenerateAugmentedSigmaPointMatrix ();

    /**
     * Run the test for the sigma point prediction step.
     */
    void testGenerateSigmaPointPrediction();

    /**
     * Run the test to predict the mean and covariance during the update step (from the predicted sigma points).
     */
    void testPredictMeanAndCovariance();

private:

    /**
     * Helper function for determining if two doubles are almost equal.
     * @param a [in] input number one.
     * @param b [in] input number two.
     * @return are the two vales within a distance of SMALL_VALUE of each other.
     */
    bool almostEqual (double a, double b) const;

    static const double SMALL_VALUE;
};