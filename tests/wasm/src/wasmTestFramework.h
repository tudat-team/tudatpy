/*    Copyright (c) 2010-2024, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    WASM test framework - common utilities and declarations.
 */

#ifndef WASM_TEST_FRAMEWORK_H
#define WASM_TEST_FRAMEWORK_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <Eigen/Core>

// Define Vector6d if not already defined
namespace Eigen {
    typedef Matrix<double, 6, 1> Vector6d;
}

// Test counters (defined in main.cpp)
extern int testsRun;
extern int testsPassed;
extern int testsFailed;

// Test assertion utilities
inline void checkClose(const std::string& testName, double actual, double expected, double tolerance = 1e-10)
{
    testsRun++;
    double diff = std::abs(actual - expected);
    bool passed = diff < tolerance;

    if (passed) {
        testsPassed++;
        std::cout << "[PASS] " << testName << std::endl;
    } else {
        testsFailed++;
        std::cout << "[FAIL] " << testName << std::endl;
        std::cout << "       Expected: " << std::setprecision(15) << expected << std::endl;
        std::cout << "       Actual:   " << std::setprecision(15) << actual << std::endl;
        std::cout << "       Diff:     " << diff << " (tolerance: " << tolerance << ")" << std::endl;
    }
}

inline void checkVectorClose(const std::string& testName, const Eigen::Vector3d& actual,
                      const Eigen::Vector3d& expected, double tolerance = 1e-10)
{
    testsRun++;
    double diff = (actual - expected).norm();
    bool passed = diff < tolerance;

    if (passed) {
        testsPassed++;
        std::cout << "[PASS] " << testName << std::endl;
    } else {
        testsFailed++;
        std::cout << "[FAIL] " << testName << std::endl;
        std::cout << "       Expected: [" << expected.transpose() << "]" << std::endl;
        std::cout << "       Actual:   [" << actual.transpose() << "]" << std::endl;
        std::cout << "       Diff norm: " << diff << " (tolerance: " << tolerance << ")" << std::endl;
    }
}

inline void checkTrue(const std::string& testName, bool condition)
{
    testsRun++;
    if (condition) {
        testsPassed++;
        std::cout << "[PASS] " << testName << std::endl;
    } else {
        testsFailed++;
        std::cout << "[FAIL] " << testName << std::endl;
    }
}

inline void checkStringEquals(const std::string& testName, const std::string& actual, const std::string& expected)
{
    testsRun++;
    if (actual == expected) {
        testsPassed++;
        std::cout << "[PASS] " << testName << std::endl;
    } else {
        testsFailed++;
        std::cout << "[FAIL] " << testName << std::endl;
        std::cout << "       Expected: \"" << expected << "\"" << std::endl;
        std::cout << "       Actual:   \"" << actual << "\"" << std::endl;
    }
}

inline void checkStringStartsWith(const std::string& testName, const std::string& actual, const std::string& prefix)
{
    testsRun++;
    bool passed = actual.substr(0, prefix.size()) == prefix;
    if (passed) {
        testsPassed++;
        std::cout << "[PASS] " << testName << std::endl;
    } else {
        testsFailed++;
        std::cout << "[FAIL] " << testName << std::endl;
        std::cout << "       Expected to start with: \"" << prefix << "\"" << std::endl;
        std::cout << "       Actual: \"" << actual << "\"" << std::endl;
    }
}

inline void checkVector6dClose(const std::string& testName, const Eigen::Vector6d& actual,
                        const Eigen::Vector6d& expected, double tolerance = 1e-10)
{
    testsRun++;
    double diff = (actual - expected).norm();
    bool passed = diff < tolerance;

    if (passed) {
        testsPassed++;
        std::cout << "[PASS] " << testName << std::endl;
    } else {
        testsFailed++;
        std::cout << "[FAIL] " << testName << std::endl;
        std::cout << "       Expected: [" << expected.transpose() << "]" << std::endl;
        std::cout << "       Actual:   [" << actual.transpose() << "]" << std::endl;
        std::cout << "       Diff norm: " << diff << " (tolerance: " << tolerance << ")" << std::endl;
    }
}

#endif // WASM_TEST_FRAMEWORK_H
