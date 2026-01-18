/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Eigen type converters for Emscripten Embind bindings.
 *    Provides JavaScript-compatible wrappers for Eigen vectors and matrices.
 */

#ifndef TUDATPY_WASM_EIGEN_H
#define TUDATPY_WASM_EIGEN_H

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <Eigen/Core>
#include <Eigen/LU>
#include <tudat/basics/basicTypedefs.h>
#include <vector>

namespace tudatpy_wasm {

// ============================================================================
// Conversion Functions: Eigen <-> JavaScript Arrays
// ============================================================================

/**
 * Convert Eigen::VectorXd to JavaScript array
 */
inline emscripten::val eigenVectorToArray(const Eigen::VectorXd& vec) {
    emscripten::val arr = emscripten::val::array();
    for (Eigen::Index i = 0; i < vec.size(); ++i) {
        arr.call<void>("push", vec(i));
    }
    return arr;
}

/**
 * Convert JavaScript array to Eigen::VectorXd
 */
inline Eigen::VectorXd arrayToEigenVector(const emscripten::val& arr) {
    unsigned length = arr["length"].as<unsigned>();
    Eigen::VectorXd vec(length);
    for (unsigned i = 0; i < length; ++i) {
        vec(i) = arr[i].as<double>();
    }
    return vec;
}

/**
 * Convert Eigen::MatrixXd to nested JavaScript array
 */
inline emscripten::val eigenMatrixToArray(const Eigen::MatrixXd& mat) {
    emscripten::val arr = emscripten::val::array();
    for (Eigen::Index i = 0; i < mat.rows(); ++i) {
        emscripten::val row = emscripten::val::array();
        for (Eigen::Index j = 0; j < mat.cols(); ++j) {
            row.call<void>("push", mat(i, j));
        }
        arr.call<void>("push", row);
    }
    return arr;
}

/**
 * Convert nested JavaScript array to Eigen::MatrixXd
 */
inline Eigen::MatrixXd arrayToEigenMatrix(const emscripten::val& arr) {
    unsigned rows = arr["length"].as<unsigned>();
    if (rows == 0) return Eigen::MatrixXd(0, 0);

    unsigned cols = arr[0]["length"].as<unsigned>();
    Eigen::MatrixXd mat(rows, cols);

    for (unsigned i = 0; i < rows; ++i) {
        for (unsigned j = 0; j < cols; ++j) {
            mat(i, j) = arr[i][j].as<double>();
        }
    }
    return mat;
}

// ============================================================================
// Fixed-Size Vector Wrappers
// ============================================================================

/**
 * Wrapper for Eigen::Vector3d with JavaScript-friendly interface
 */
struct Vector3dWrapper {
    Eigen::Vector3d data;

    Vector3dWrapper() : data(Eigen::Vector3d::Zero()) {}
    Vector3dWrapper(double x, double y, double z) : data(x, y, z) {}
    explicit Vector3dWrapper(const Eigen::Vector3d& v) : data(v) {}

    double get(int i) const { return data(i); }
    void set(int i, double val) { data(i) = val; }

    double x() const { return data(0); }
    double y() const { return data(1); }
    double z() const { return data(2); }

    void setX(double val) { data(0) = val; }
    void setY(double val) { data(1) = val; }
    void setZ(double val) { data(2) = val; }

    double norm() const { return data.norm(); }
    Vector3dWrapper normalized() const { return Vector3dWrapper(data.normalized()); }

    emscripten::val toArray() const {
        emscripten::val arr = emscripten::val::array();
        arr.call<void>("push", data(0));
        arr.call<void>("push", data(1));
        arr.call<void>("push", data(2));
        return arr;
    }

    static Vector3dWrapper fromArray(const emscripten::val& arr) {
        return Vector3dWrapper(
            arr[0].as<double>(),
            arr[1].as<double>(),
            arr[2].as<double>()
        );
    }

    const Eigen::Vector3d& eigen() const { return data; }
};

/**
 * Wrapper for Eigen::Vector6d with JavaScript-friendly interface
 */
struct Vector6dWrapper {
    Eigen::Vector6d data;

    Vector6dWrapper() : data(Eigen::Vector6d::Zero()) {}
    Vector6dWrapper(double a, double b, double c, double d, double e, double f) {
        data << a, b, c, d, e, f;
    }
    explicit Vector6dWrapper(const Eigen::Vector6d& v) : data(v) {}

    double get(int i) const { return data(i); }
    void set(int i, double val) { data(i) = val; }

    int size() const { return 6; }
    double norm() const { return data.norm(); }

    Vector3dWrapper position() const { return Vector3dWrapper(data.head<3>()); }
    Vector3dWrapper velocity() const { return Vector3dWrapper(data.tail<3>()); }

    emscripten::val toArray() const {
        emscripten::val arr = emscripten::val::array();
        for (int i = 0; i < 6; ++i) {
            arr.call<void>("push", data(i));
        }
        return arr;
    }

    static Vector6dWrapper fromArray(const emscripten::val& arr) {
        Vector6dWrapper w;
        for (int i = 0; i < 6; ++i) {
            w.data(i) = arr[i].as<double>();
        }
        return w;
    }

    const Eigen::Vector6d& eigen() const { return data; }
};

/**
 * Wrapper for Eigen::Vector7d (quaternion + angular velocity)
 */
struct Vector7dWrapper {
    Eigen::Matrix<double, 7, 1> data;

    Vector7dWrapper() : data(Eigen::Matrix<double, 7, 1>::Zero()) {}
    explicit Vector7dWrapper(const Eigen::Matrix<double, 7, 1>& v) : data(v) {}

    double get(int i) const { return data(i); }
    void set(int i, double val) { data(i) = val; }
    int size() const { return 7; }

    emscripten::val toArray() const {
        emscripten::val arr = emscripten::val::array();
        for (int i = 0; i < 7; ++i) {
            arr.call<void>("push", data(i));
        }
        return arr;
    }

    static Vector7dWrapper fromArray(const emscripten::val& arr) {
        Vector7dWrapper w;
        for (int i = 0; i < 7; ++i) {
            w.data(i) = arr[i].as<double>();
        }
        return w;
    }
};

// ============================================================================
// Dynamic Vector Wrapper
// ============================================================================

/**
 * Wrapper for Eigen::VectorXd with JavaScript-friendly interface
 */
struct VectorXdWrapper {
    Eigen::VectorXd data;

    VectorXdWrapper() : data() {}
    explicit VectorXdWrapper(int n) : data(Eigen::VectorXd::Zero(n)) {}
    explicit VectorXdWrapper(const Eigen::VectorXd& v) : data(v) {}

    double get(int i) const { return data(i); }
    void set(int i, double val) { data(i) = val; }
    int size() const { return static_cast<int>(data.size()); }
    void resize(int n) { data.resize(n); }
    double norm() const { return data.norm(); }

    emscripten::val toArray() const { return eigenVectorToArray(data); }

    static VectorXdWrapper fromArray(const emscripten::val& arr) {
        return VectorXdWrapper(arrayToEigenVector(arr));
    }

    const Eigen::VectorXd& eigen() const { return data; }
};

// ============================================================================
// Matrix Wrappers
// ============================================================================

/**
 * Wrapper for Eigen::Matrix3d
 */
struct Matrix3dWrapper {
    Eigen::Matrix3d data;

    Matrix3dWrapper() : data(Eigen::Matrix3d::Zero()) {}
    explicit Matrix3dWrapper(const Eigen::Matrix3d& m) : data(m) {}

    static Matrix3dWrapper identity() { return Matrix3dWrapper(Eigen::Matrix3d::Identity()); }

    double get(int r, int c) const { return data(r, c); }
    void set(int r, int c, double val) { data(r, c) = val; }

    int rows() const { return 3; }
    int cols() const { return 3; }

    Matrix3dWrapper transpose() const { return Matrix3dWrapper(data.transpose()); }
    Matrix3dWrapper inverse() const { return Matrix3dWrapper(data.inverse()); }
    double determinant() const { return data.determinant(); }

    Vector3dWrapper multiply(const Vector3dWrapper& v) const {
        return Vector3dWrapper(data * v.data);
    }

    emscripten::val toArray() const {
        emscripten::val arr = emscripten::val::array();
        for (int i = 0; i < 3; ++i) {
            emscripten::val row = emscripten::val::array();
            for (int j = 0; j < 3; ++j) {
                row.call<void>("push", data(i, j));
            }
            arr.call<void>("push", row);
        }
        return arr;
    }

    const Eigen::Matrix3d& eigen() const { return data; }
};

/**
 * Wrapper for Eigen::MatrixXd
 */
struct MatrixXdWrapper {
    Eigen::MatrixXd data;

    MatrixXdWrapper() : data() {}
    MatrixXdWrapper(int r, int c) : data(Eigen::MatrixXd::Zero(r, c)) {}
    explicit MatrixXdWrapper(const Eigen::MatrixXd& m) : data(m) {}

    double get(int r, int c) const { return data(r, c); }
    void set(int r, int c, double val) { data(r, c) = val; }

    int rows() const { return static_cast<int>(data.rows()); }
    int cols() const { return static_cast<int>(data.cols()); }
    void resize(int r, int c) { data.resize(r, c); }

    MatrixXdWrapper transpose() const { return MatrixXdWrapper(data.transpose()); }

    emscripten::val toArray() const { return eigenMatrixToArray(data); }

    static MatrixXdWrapper fromArray(const emscripten::val& arr) {
        return MatrixXdWrapper(arrayToEigenMatrix(arr));
    }

    const Eigen::MatrixXd& eigen() const { return data; }
};

} // namespace tudatpy_wasm

// Embind registrations are in eigen_wasm.cpp (not in header to avoid duplicate symbols)

#endif // __EMSCRIPTEN__

#endif // TUDATPY_WASM_EIGEN_H
