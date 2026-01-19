/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include "../../wasm_module.h"
#include "../../eigen_wasm.h"

#include <tudat/astro/basic_astro/polyhedronFuntions.h>

namespace tba = tudat::basic_astrodynamics;

WASM_MODULE_PATH("astro_polyhedron_utilities")

// ============================================================================
// Wrapper functions for polyhedron utilities
// These convert between WASM wrapper types and native Eigen types
// ============================================================================

// Helper function to convert MatrixXdWrapper to Eigen::MatrixXd
Eigen::MatrixXd matrixXdWrapperToEigen(const tudatpy_wasm::MatrixXdWrapper& wrapper) {
    return wrapper.data;
}

// Helper function to convert MatrixXi-like data (stored as MatrixXd in JS) to Eigen::MatrixXi
// JavaScript doesn't have native integer matrices, so we receive doubles and convert
Eigen::MatrixXi matrixXdToMatrixXi(const tudatpy_wasm::MatrixXdWrapper& wrapper) {
    return wrapper.data.cast<int>();
}

// Wrapper: Compute polyhedron surface area
double polyhedronSurfaceAreaWrapper(
    const tudatpy_wasm::MatrixXdWrapper& verticesCoordinates,
    const tudatpy_wasm::MatrixXdWrapper& verticesDefiningEachFacet)
{
    return tba::computePolyhedronSurfaceArea(
        matrixXdWrapperToEigen(verticesCoordinates),
        matrixXdToMatrixXi(verticesDefiningEachFacet));
}

// Wrapper: Compute polyhedron volume
double polyhedronVolumeWrapper(
    const tudatpy_wasm::MatrixXdWrapper& verticesCoordinates,
    const tudatpy_wasm::MatrixXdWrapper& verticesDefiningEachFacet)
{
    return tba::computePolyhedronVolume(
        matrixXdWrapperToEigen(verticesCoordinates),
        matrixXdToMatrixXi(verticesDefiningEachFacet));
}

// Wrapper: Compute polyhedron mean radius
double polyhedronMeanRadiusWrapper(
    const tudatpy_wasm::MatrixXdWrapper& verticesCoordinates,
    const tudatpy_wasm::MatrixXdWrapper& verticesDefiningEachFacet)
{
    return tba::computePolyhedronMeanRadius(
        matrixXdWrapperToEigen(verticesCoordinates),
        matrixXdToMatrixXi(verticesDefiningEachFacet));
}

// Wrapper: Compute polyhedron centroid position
tudatpy_wasm::Vector3dWrapper polyhedronCentroidWrapper(
    const tudatpy_wasm::MatrixXdWrapper& verticesCoordinates,
    const tudatpy_wasm::MatrixXdWrapper& verticesDefiningEachFacet)
{
    Eigen::Vector3d centroid = tba::computePolyhedronCentroidPosition(
        matrixXdWrapperToEigen(verticesCoordinates),
        matrixXdToMatrixXi(verticesDefiningEachFacet));
    return tudatpy_wasm::Vector3dWrapper(centroid);
}

// Wrapper: Modify polyhedron centroid position
tudatpy_wasm::MatrixXdWrapper polyhedronModifyCentroidWrapper(
    const tudatpy_wasm::MatrixXdWrapper& verticesCoordinates,
    const tudatpy_wasm::MatrixXdWrapper& verticesDefiningEachFacet,
    const tudatpy_wasm::Vector3dWrapper& desiredCentroid)
{
    Eigen::MatrixXd result = tba::modifyPolyhedronCentroidPosition(
        matrixXdWrapperToEigen(verticesCoordinates),
        matrixXdToMatrixXi(verticesDefiningEachFacet),
        desiredCentroid.data);
    return tudatpy_wasm::MatrixXdWrapper(result);
}

// Wrapper: Compute polyhedron inertia tensor from density
tudatpy_wasm::Matrix3dWrapper polyhedronInertiaTensorFromDensityWrapper(
    const tudatpy_wasm::MatrixXdWrapper& verticesCoordinates,
    const tudatpy_wasm::MatrixXdWrapper& verticesDefiningEachFacet,
    double density)
{
    Eigen::Matrix3d tensor = tba::computePolyhedronInertiaTensor(
        matrixXdWrapperToEigen(verticesCoordinates),
        matrixXdToMatrixXi(verticesDefiningEachFacet),
        density);
    return tudatpy_wasm::Matrix3dWrapper(tensor);
}

// Wrapper: Compute polyhedron inertia tensor from gravitational parameter
tudatpy_wasm::Matrix3dWrapper polyhedronInertiaTensorFromGravParamWrapper(
    const tudatpy_wasm::MatrixXdWrapper& verticesCoordinates,
    const tudatpy_wasm::MatrixXdWrapper& verticesDefiningEachFacet,
    double gravitationalParameter,
    double gravitationalConstant)
{
    Eigen::Matrix3d tensor = tba::computePolyhedronInertiaTensor(
        matrixXdWrapperToEigen(verticesCoordinates),
        matrixXdToMatrixXi(verticesDefiningEachFacet),
        gravitationalParameter,
        gravitationalConstant);
    return tudatpy_wasm::Matrix3dWrapper(tensor);
}

EMSCRIPTEN_BINDINGS(tudatpy_astro_polyhedron_utilities) {
    using namespace emscripten;

    // ========================================================================
    // Polyhedron utility functions
    // ========================================================================

    // Surface area computation
    function("astro_polyhedron_utilities_surface_area",
        &polyhedronSurfaceAreaWrapper);

    // Volume computation
    function("astro_polyhedron_utilities_volume",
        &polyhedronVolumeWrapper);

    // Mean radius computation
    function("astro_polyhedron_utilities_mean_radius",
        &polyhedronMeanRadiusWrapper);

    // Centroid position computation
    function("astro_polyhedron_utilities_centroid",
        &polyhedronCentroidWrapper);

    // Modify centroid position
    function("astro_polyhedron_utilities_modify_centroid",
        &polyhedronModifyCentroidWrapper);

    // Inertia tensor from density
    function("astro_polyhedron_utilities_inertia_tensor_from_density",
        &polyhedronInertiaTensorFromDensityWrapper);

    // Inertia tensor from gravitational parameter
    function("astro_polyhedron_utilities_inertia_tensor_from_gravitational_parameter",
        &polyhedronInertiaTensorFromGravParamWrapper);
}

#endif
