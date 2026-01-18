/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Embind registration for Eigen type wrappers.
 *    This is separate from the header to avoid duplicate symbol errors.
 */

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include "eigen_wasm.h"

// ============================================================================
// Embind Registration for Eigen Types
// ============================================================================

EMSCRIPTEN_BINDINGS(eigen_types) {
    using namespace emscripten;
    using namespace tudatpy_wasm;

    // Vector3d
    class_<Vector3dWrapper>("Vector3d")
        .constructor<>()
        .constructor<double, double, double>()
        .function("get", &Vector3dWrapper::get)
        .function("set", &Vector3dWrapper::set)
        .function("size", &Vector3dWrapper::size)
        .property("x", &Vector3dWrapper::x, &Vector3dWrapper::setX)
        .property("y", &Vector3dWrapper::y, &Vector3dWrapper::setY)
        .property("z", &Vector3dWrapper::z, &Vector3dWrapper::setZ)
        .function("norm", &Vector3dWrapper::norm)
        .function("normalized", &Vector3dWrapper::normalized)
        .function("toArray", &Vector3dWrapper::toArray)
        .class_function("fromArray", &Vector3dWrapper::fromArray);

    // Vector6d
    class_<Vector6dWrapper>("Vector6d")
        .constructor<>()
        .constructor<double, double, double, double, double, double>()
        .function("get", &Vector6dWrapper::get)
        .function("set", &Vector6dWrapper::set)
        .function("size", &Vector6dWrapper::size)
        .function("norm", &Vector6dWrapper::norm)
        .function("position", &Vector6dWrapper::position)
        .function("velocity", &Vector6dWrapper::velocity)
        .function("toArray", &Vector6dWrapper::toArray)
        .class_function("fromArray", &Vector6dWrapper::fromArray);

    // Vector7d
    class_<Vector7dWrapper>("Vector7d")
        .constructor<>()
        .function("get", &Vector7dWrapper::get)
        .function("set", &Vector7dWrapper::set)
        .function("size", &Vector7dWrapper::size)
        .function("toArray", &Vector7dWrapper::toArray)
        .class_function("fromArray", &Vector7dWrapper::fromArray);

    // VectorXd
    class_<VectorXdWrapper>("VectorXd")
        .constructor<>()
        .constructor<int>()
        .function("get", &VectorXdWrapper::get)
        .function("set", &VectorXdWrapper::set)
        .function("size", &VectorXdWrapper::size)
        .function("resize", &VectorXdWrapper::resize)
        .function("norm", &VectorXdWrapper::norm)
        .function("toArray", &VectorXdWrapper::toArray)
        .class_function("fromArray", &VectorXdWrapper::fromArray);

    // Matrix3d
    class_<Matrix3dWrapper>("Matrix3d")
        .constructor<>()
        .class_function("identity", &Matrix3dWrapper::identity)
        .function("get", &Matrix3dWrapper::get)
        .function("set", &Matrix3dWrapper::set)
        .function("rows", &Matrix3dWrapper::rows)
        .function("cols", &Matrix3dWrapper::cols)
        .function("transpose", &Matrix3dWrapper::transpose)
        .function("inverse", &Matrix3dWrapper::inverse)
        .function("determinant", &Matrix3dWrapper::determinant)
        .function("multiply", &Matrix3dWrapper::multiply)
        .function("toArray", &Matrix3dWrapper::toArray);

    // MatrixXd
    class_<MatrixXdWrapper>("MatrixXd")
        .constructor<>()
        .constructor<int, int>()
        .function("get", &MatrixXdWrapper::get)
        .function("set", &MatrixXdWrapper::set)
        .function("rows", &MatrixXdWrapper::rows)
        .function("cols", &MatrixXdWrapper::cols)
        .function("resize", &MatrixXdWrapper::resize)
        .function("transpose", &MatrixXdWrapper::transpose)
        .function("toArray", &MatrixXdWrapper::toArray)
        .class_function("fromArray", &MatrixXdWrapper::fromArray);
}

#endif // __EMSCRIPTEN__
