/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Embind registration for STL container types.
 *    This is separate from the header to avoid duplicate symbol errors.
 */

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include <vector>
#include <map>
#include <string>
#include <memory>
#include <Eigen/Core>

#include "eigen_wasm.h"

// ============================================================================
// STL Container Registrations
// ============================================================================

EMSCRIPTEN_BINDINGS(stl_containers) {
    using namespace emscripten;
    using namespace tudatpy_wasm;

    // ========================================================================
    // Vectors of Primitives
    // ========================================================================

    register_vector<double>("VectorDouble");
    register_vector<int>("VectorInt");
    register_vector<unsigned int>("VectorUInt");
    register_vector<std::string>("VectorString");
    // Note: register_vector<bool> is not supported due to std::vector<bool> specialization

    // ========================================================================
    // Vectors of Vectors (Nested)
    // ========================================================================

    register_vector<std::vector<double>>("VectorVectorDouble");

    // ========================================================================
    // Vectors of Eigen Wrappers
    // ========================================================================

    register_vector<Vector3dWrapper>("VectorVector3d");
    register_vector<Vector6dWrapper>("VectorVector6d");
    register_vector<VectorXdWrapper>("VectorVectorXd");
    register_vector<Matrix3dWrapper>("VectorMatrix3d");
    register_vector<MatrixXdWrapper>("VectorMatrixXd");

    // ========================================================================
    // Maps with String Keys
    // ========================================================================

    register_map<std::string, double>("MapStringDouble");
    register_map<std::string, int>("MapStringInt");
    register_map<std::string, std::string>("MapStringString");
    register_map<std::string, bool>("MapStringBool");

    // ========================================================================
    // Maps with Double Keys (Time-Keyed Data - Critical for Propagation)
    // ========================================================================

    register_map<double, double>("MapDoubleDouble");
    register_map<double, Vector3dWrapper>("MapDoubleVector3d");
    register_map<double, Vector6dWrapper>("MapDoubleVector6d");
    register_map<double, VectorXdWrapper>("MapDoubleVectorXd");
    register_map<double, MatrixXdWrapper>("MapDoubleMatrixXd");

    // ========================================================================
    // Maps with String Keys to Eigen Types
    // ========================================================================

    register_map<std::string, Vector3dWrapper>("MapStringVector3d");
    register_map<std::string, Vector6dWrapper>("MapStringVector6d");
    register_map<std::string, VectorXdWrapper>("MapStringVectorXd");

    // ========================================================================
    // Maps with String Keys to Vectors
    // ========================================================================

    register_map<std::string, std::vector<double>>("MapStringVectorDouble");
    register_map<std::string, std::vector<std::string>>("MapStringVectorString");

    // ========================================================================
    // Pairs (for map iteration)
    // ========================================================================

    value_object<std::pair<double, double>>("PairDoubleDouble")
        .field("first", &std::pair<double, double>::first)
        .field("second", &std::pair<double, double>::second);

    value_object<std::pair<std::string, double>>("PairStringDouble")
        .field("first", &std::pair<std::string, double>::first)
        .field("second", &std::pair<std::string, double>::second);

    value_object<std::pair<std::string, std::string>>("PairStringString")
        .field("first", &std::pair<std::string, std::string>::first)
        .field("second", &std::pair<std::string, std::string>::second);

    // ========================================================================
    // Pairs for Covariance Propagation Results
    // ========================================================================

    // Pair for split output covariance propagation (times + matrices)
    value_object<std::pair<std::vector<double>, std::vector<MatrixXdWrapper>>>(
        "PairVectorDoubleVectorMatrixXd")
        .field("times", &std::pair<std::vector<double>, std::vector<MatrixXdWrapper>>::first)
        .field("matrices", &std::pair<std::vector<double>, std::vector<MatrixXdWrapper>>::second);

    // Pair for split output formal errors propagation (times + vectors)
    value_object<std::pair<std::vector<double>, std::vector<VectorXdWrapper>>>(
        "PairVectorDoubleVectorVectorXd")
        .field("times", &std::pair<std::vector<double>, std::vector<VectorXdWrapper>>::first)
        .field("vectors", &std::pair<std::vector<double>, std::vector<VectorXdWrapper>>::second);
}

#endif // __EMSCRIPTEN__
