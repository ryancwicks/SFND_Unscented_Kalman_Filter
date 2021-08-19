/**
 * main.cpp
 * 
 * Main entry point for the test code.
 * 
 * Ryan Wicks, 17 August 2021
 */

#include "UKFTest.h"
#include <iostream>

int main (int argc, char * argv[]) {

    UKFTest test;

    test.testGenerateSigmaPointMatrix();
    test.testGenerateAugmentedSigmaPointMatrix();
    test.testGenerateSigmaPointPrediction();
    return 0;
}

