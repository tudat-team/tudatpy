/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Element conversion module bindings for WASM.
 *    Mirrors: src/tudatpy/astro/element_conversion/expose_element_conversion.cpp
 *
 *    Provides orbital element conversions between different representations:
 *    - Cartesian <-> Keplerian
 *    - Spherical <-> Cartesian
 *    - Modified Equinoctial Elements (MEE)
 *    - Unified State Model (USM)
 *    - And more
 */

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include "expose_element_conversion_wasm.h"
#include "../../wasm_module.h"
#include "../../eigen_wasm.h"

// Tudat headers
#include <tudat/astro/basic_astro/attitudeElementConversions.h>
#include <tudat/astro/basic_astro/stateRepresentationConversions.h>
#include <tudat/astro/conversions.h>
#include <tudat/astro/ephemerides/rotationalEphemeris.h>
#include <tudat/math/basic.h>
#include <tudat/interface/spice/spiceInterface.h>
#include <tudat/astro/basic_astro/missionGeometry.h>

namespace toec = tudat::orbital_element_conversions;
namespace tcc = tudat::coordinate_conversions;
namespace tla = tudat::linear_algebra;
namespace te = tudat::ephemerides;
namespace tba = tudat::basic_astrodynamics;
namespace tmg = tudat::mission_geometry;
namespace tsi = tudat::spice_interface;

// Register this module path
WASM_MODULE_PATH("astro_element_conversion")

// ============================================================================
// Wrapper Functions for Eigen Vector Conversions
// These adapt between JavaScript-friendly wrappers and native Eigen types
// ============================================================================

namespace {

using namespace tudatpy_wasm;

// Keplerian to Cartesian conversion wrapper
Vector6dWrapper keplerianToCartesian(const Vector6dWrapper& keplerian, double gravitationalParameter) {
    return Vector6dWrapper(toec::convertKeplerianToCartesianElements<double>(
        keplerian.eigen(), gravitationalParameter));
}

// Cartesian to Keplerian conversion wrapper
Vector6dWrapper cartesianToKeplerian(const Vector6dWrapper& cartesian, double gravitationalParameter) {
    return Vector6dWrapper(toec::convertCartesianToKeplerianElements<double>(
        cartesian.eigen(), gravitationalParameter));
}

// Spherical to Cartesian state conversion wrapper
Vector6dWrapper sphericalToCartesianState(const Vector6dWrapper& spherical) {
    return Vector6dWrapper(tcc::convertSphericalToCartesianState(spherical.eigen()));
}

// Cartesian to Spherical state conversion wrapper
Vector6dWrapper cartesianToSphericalState(const Vector6dWrapper& cartesian) {
    return Vector6dWrapper(tcc::convertCartesianToSphericalState(cartesian.eigen()));
}

// MEE to Cartesian wrapper
Vector6dWrapper meeToCartesian(const Vector6dWrapper& mee, double gravitationalParameter, bool flipSingularity) {
    return Vector6dWrapper(toec::convertModifiedEquinoctialToCartesianElements<double>(
        mee.eigen(), gravitationalParameter, flipSingularity));
}

// Cartesian to MEE wrapper
Vector6dWrapper cartesianToMee(const Vector6dWrapper& cartesian, double gravitationalParameter, bool flipSingularity) {
    return Vector6dWrapper(toec::convertCartesianToModifiedEquinoctialElements<double>(
        cartesian.eigen(), gravitationalParameter, flipSingularity));
}

// Keplerian to MEE wrapper
Vector6dWrapper keplerianToMee(const Vector6dWrapper& keplerian, bool flipSingularity) {
    return Vector6dWrapper(toec::convertKeplerianToModifiedEquinoctialElements<double>(
        keplerian.eigen(), flipSingularity));
}

// MEE to Keplerian wrapper
Vector6dWrapper meeToKeplerian(const Vector6dWrapper& mee, bool flipSingularity) {
    return Vector6dWrapper(toec::convertModifiedEquinoctialToKeplerianElements<double>(
        mee.eigen(), flipSingularity));
}

// Geodetic coordinate conversions
Vector3dWrapper cartesianToGeodetic(const Vector3dWrapper& cartesian,
                                     double equatorialRadius, double flattening, double tolerance) {
    return Vector3dWrapper(tcc::convertCartesianToGeodeticCoordinates(
        cartesian.eigen(), equatorialRadius, flattening, tolerance));
}

Vector3dWrapper geodeticToCartesian(const Vector3dWrapper& geodetic,
                                     double equatorialRadius, double flattening) {
    return Vector3dWrapper(tcc::convertGeodeticToCartesianCoordinates(
        geodetic.eigen(), equatorialRadius, flattening));
}

// Spherical position conversions
Vector3dWrapper cartesianToSphericalPosition(const Vector3dWrapper& cartesian) {
    return Vector3dWrapper(tcc::convertCartesianToSpherical(cartesian.eigen()));
}

Vector3dWrapper sphericalToCartesianPosition(const Vector3dWrapper& spherical) {
    return Vector3dWrapper(tcc::convertSphericalToCartesian(spherical.eigen()));
}

// True anomaly conversions
double meanToTrueAnomaly(double eccentricity, double meanAnomaly) {
    return toec::convertMeanAnomalyToTrueAnomaly<double>(eccentricity, meanAnomaly);
}

double trueToMeanAnomaly(double eccentricity, double trueAnomaly) {
    return toec::convertTrueAnomalyToMeanAnomaly<double>(eccentricity, trueAnomaly);
}

double meanToEccentricAnomaly(double eccentricity, double meanAnomaly) {
    return toec::convertMeanAnomalyToEccentricAnomaly<double>(eccentricity, meanAnomaly);
}

double eccentricToMeanAnomaly(double eccentricity, double eccentricAnomaly) {
    return toec::convertEccentricAnomalyToMeanAnomaly<double>(eccentricity, eccentricAnomaly);
}

double eccentricToTrueAnomaly(double eccentricAnomaly, double eccentricity) {
    return toec::convertEccentricAnomalyToTrueAnomaly<double>(eccentricAnomaly, eccentricity);
}

double trueToEccentricAnomaly(double trueAnomaly, double eccentricity) {
    return toec::convertTrueAnomalyToEccentricAnomaly<double>(trueAnomaly, eccentricity);
}

// Hyperbolic anomaly conversions
double meanToHyperbolicEccentricAnomaly(double eccentricity, double hyperbolicMeanAnomaly) {
    return toec::convertMeanAnomalyToHyperbolicEccentricAnomaly<double>(eccentricity, hyperbolicMeanAnomaly);
}

double hyperbolicEccentricToMeanAnomaly(double eccentricity, double hyperbolicEccentricAnomaly) {
    return toec::convertHyperbolicEccentricAnomalyToMeanAnomaly<double>(eccentricity, hyperbolicEccentricAnomaly);
}

double hyperbolicEccentricToTrueAnomaly(double hyperbolicEccentricAnomaly, double eccentricity) {
    return toec::convertHyperbolicEccentricAnomalyToTrueAnomaly<double>(hyperbolicEccentricAnomaly, eccentricity);
}

double trueToHyperbolicEccentricAnomaly(double trueAnomaly, double eccentricity) {
    return toec::convertTrueAnomalyToHyperbolicEccentricAnomaly<double>(trueAnomaly, eccentricity);
}

// USM element conversions (all USM formats return 7-element vectors)
Vector7dWrapper cartesianToUsmEm(const Vector6dWrapper& cartesian, double gravitationalParameter) {
    Eigen::Matrix<double, 7, 1> result = toec::convertCartesianToUnifiedStateModelExponentialMapElements(
        cartesian.data, gravitationalParameter);
    return Vector7dWrapper(result);
}

Vector7dWrapper cartesianToUsm7(const Vector6dWrapper& cartesian, double gravitationalParameter) {
    Eigen::Matrix<double, 7, 1> result = toec::convertCartesianToUnifiedStateModelQuaternionsElements(
        cartesian.data, gravitationalParameter);
    return Vector7dWrapper(result);
}

Vector7dWrapper cartesianToUsm6(const Vector6dWrapper& cartesian, double gravitationalParameter) {
    Eigen::Matrix<double, 7, 1> result = toec::convertCartesianToUnifiedStateModelModifiedRodriguesParameterElements(
        cartesian.data, gravitationalParameter);
    return Vector7dWrapper(result);
}

Vector6dWrapper usmEmToCartesian(const Vector7dWrapper& usm, double gravitationalParameter) {
    return Vector6dWrapper(toec::convertUnifiedStateModelExponentialMapToCartesianElements(
        usm.data, gravitationalParameter));
}

Vector6dWrapper usm7ToCartesian(const Vector7dWrapper& usm, double gravitationalParameter, bool normalizeQuaternion = true) {
    return Vector6dWrapper(toec::convertUnifiedStateModelQuaternionsToCartesianElements(
        usm.data, gravitationalParameter, normalizeQuaternion));
}

Vector6dWrapper usm6ToCartesian(const Vector7dWrapper& usm, double gravitationalParameter) {
    return Vector6dWrapper(toec::convertUnifiedStateModelModifiedRodriguesParametersToCartesianElements(
        usm.data, gravitationalParameter));
}

// Mean motion / elapsed time conversions
double elapsedTimeToDeltaMeanAnomaly(double elapsedTime, double gravitationalParameter, double semiMajorAxis) {
    return toec::convertElapsedTimeToMeanAnomalyChange<double>(elapsedTime, gravitationalParameter, semiMajorAxis);
}

double deltaMeanAnomalyToElapsedTime(double meanAnomalyChange, double gravitationalParameter, double semiMajorAxis) {
    return toec::convertMeanAnomalyChangeToElapsedTime<double>(meanAnomalyChange, gravitationalParameter, semiMajorAxis);
}

double meanMotionToSemiMajorAxis(double meanMotion, double gravitationalParameter) {
    return toec::convertEllipticalMeanMotionToSemiMajorAxis<double>(meanMotion, gravitationalParameter);
}

double semiMajorAxisToMeanMotion(double semiMajorAxis, double gravitationalParameter) {
    return toec::convertSemiMajorAxisToEllipticalMeanMotion<double>(semiMajorAxis, gravitationalParameter);
}

// Quaternion conversions
Matrix3dWrapper quaternionEntriesToRotationMatrix(const VectorXdWrapper& quaternionEntries) {
    return Matrix3dWrapper(tla::convertVectorQuaternionToMatrixFormat(quaternionEntries.data));
}

VectorXdWrapper rotationMatrixToQuaternionEntries(const Matrix3dWrapper& rotationMatrix) {
    return VectorXdWrapper(tla::convertMatrixToVectorQuaternionFormat(rotationMatrix.data));
}

// Quaternion <-> MRP conversions (4-element inputs/outputs)
VectorXdWrapper quaternionToModifiedRodriguesParameters(const VectorXdWrapper& quaternionEntries) {
    // Takes 4-element quaternion, returns 4-element MRP (3 MRP elements + shadow flag)
    Eigen::Vector4d quat = quaternionEntries.data.head<4>();
    Eigen::Vector4d result = toec::convertQuaternionsToModifiedRodriguesParameterElements(quat);
    return VectorXdWrapper(result);
}

VectorXdWrapper modifiedRodriguesParametersToQuaternion(const VectorXdWrapper& mrp) {
    // Takes 4-element MRP (3 MRP elements + shadow flag), returns 4-element quaternion
    Eigen::Vector4d mrpVec = mrp.data.head<4>();
    Eigen::Vector4d result = toec::convertModifiedRodriguesParametersToQuaternionElements(mrpVec);
    return VectorXdWrapper(result);
}

// Quaternion <-> Exponential map conversions (4-element inputs/outputs)
VectorXdWrapper quaternionToExponentialMap(const VectorXdWrapper& quaternionEntries) {
    // Takes 4-element quaternion, returns 4-element exponential map
    Eigen::Vector4d quat = quaternionEntries.data.head<4>();
    Eigen::Vector4d result = toec::convertQuaternionsToExponentialMapElements(quat);
    return VectorXdWrapper(result);
}

VectorXdWrapper exponentialMapToQuaternion(const VectorXdWrapper& expMap) {
    // Takes 4-element exponential map, returns 4-element quaternion
    Eigen::Vector4d expMapVec = expMap.data.head<4>();
    Eigen::Vector4d result = toec::convertExponentialMapToQuaternionElements(expMapVec);
    return VectorXdWrapper(result);
}

// TEME/J2000 conversions
Matrix3dWrapper temeToJ2000(double epoch) {
    return Matrix3dWrapper(te::getRotationMatrixFromTemeToJ2000(epoch));
}

Matrix3dWrapper j2000ToTeme(double epoch) {
    // Inverse is transpose for rotation matrices
    return Matrix3dWrapper(te::getRotationMatrixFromTemeToJ2000(epoch).transpose());
}

// J2000/ECLIPJ2000 conversions
Matrix3dWrapper j2000ToEclipJ2000() {
    return Matrix3dWrapper(tsi::getRotationFromJ2000ToEclipJ2000());
}

Matrix3dWrapper eclipJ2000ToJ2000() {
    return Matrix3dWrapper(tsi::getRotationFromEclipJ2000ToJ2000());
}

// Keplerian to Cartesian elementwise
Vector6dWrapper keplerianToCartesianElementwise(double semiMajorAxis, double eccentricity, double inclination,
                                                 double argumentOfPeriapsis, double longitudeOfAscendingNode,
                                                 double trueAnomaly, double gravitationalParameter) {
    return Vector6dWrapper(toec::convertKeplerianToCartesianElements<double>(
        semiMajorAxis, eccentricity, inclination, argumentOfPeriapsis, longitudeOfAscendingNode,
        trueAnomaly, gravitationalParameter));
}

// Spherical orbital to Cartesian elementwise
Vector6dWrapper sphericalToCartesianElementwise(double radialDistance, double latitude, double longitude,
                                                 double speed, double flightPathAngle, double headingAngle) {
    return Vector6dWrapper(toec::convertSphericalOrbitalToCartesianState<double>(
        radialDistance, latitude, longitude, speed, flightPathAngle, headingAngle));
}

// Spherical orbital <-> Cartesian (full state)
Vector6dWrapper sphericalOrbitalToCartesian(const Vector6dWrapper& spherical) {
    return Vector6dWrapper(toec::convertSphericalOrbitalToCartesianState<double>(spherical.data));
}

Vector6dWrapper cartesianToSphericalOrbital(const Vector6dWrapper& cartesian) {
    return Vector6dWrapper(toec::convertCartesianToSphericalOrbitalState(cartesian.data));
}

// Flip MEE singularity helper
bool flipMeeSingularity(const Vector6dWrapper& keplerianElements) {
    return tmg::isOrbitRetrograde(keplerianElements.data);
}

// Keplerian to MEE (auto singularity)
Vector6dWrapper keplerianToMeeAuto(const Vector6dWrapper& keplerian) {
    return Vector6dWrapper(toec::convertKeplerianToModifiedEquinoctialElements<double>(keplerian.data));
}

// Cartesian to MEE (auto singularity)
Vector6dWrapper cartesianToMeeAuto(const Vector6dWrapper& cartesian, double gravitationalParameter) {
    return Vector6dWrapper(toec::convertCartesianToModifiedEquinoctialElements<double>(
        cartesian.data, gravitationalParameter));
}

} // anonymous namespace

// ============================================================================
// Embind Registration
// ============================================================================

EMSCRIPTEN_BINDINGS(tudatpy_astro_element_conversion) {
    using namespace emscripten;

    // ========================================================================
    // Enumerations
    // ========================================================================

    enum_<toec::KeplerianElementIndices>("astro_element_conversion_KeplerianElementIndices")
        .value("semi_major_axis_index", toec::KeplerianElementIndices::semiMajorAxisIndex)
        .value("semi_latus_rectum_index", toec::KeplerianElementIndices::semiLatusRectumIndex)
        .value("eccentricity_index", toec::KeplerianElementIndices::eccentricityIndex)
        .value("inclination_index", toec::KeplerianElementIndices::inclinationIndex)
        .value("argument_of_periapsis_index", toec::KeplerianElementIndices::argumentOfPeriapsisIndex)
        .value("longitude_of_ascending_node_index", toec::KeplerianElementIndices::longitudeOfAscendingNodeIndex)
        .value("true_anomaly_index", toec::KeplerianElementIndices::trueAnomalyIndex);

    enum_<toec::SphericalOrbitalStateElementIndices>("astro_element_conversion_SphericalOrbitalStateElementIndices")
        .value("radius_index", toec::SphericalOrbitalStateElementIndices::radiusIndex)
        .value("latitude_index", toec::SphericalOrbitalStateElementIndices::latitudeIndex)
        .value("longitude_index", toec::SphericalOrbitalStateElementIndices::longitudeIndex)
        .value("speed_index", toec::SphericalOrbitalStateElementIndices::speedIndex)
        .value("flight_path_index", toec::SphericalOrbitalStateElementIndices::flightPathIndex)
        .value("heading_angle_index", toec::SphericalOrbitalStateElementIndices::headingAngleIndex);

    enum_<tcc::PositionElementTypes>("astro_element_conversion_PositionElementTypes")
        .value("cartesian_position_type", tcc::PositionElementTypes::cartesian_position)
        .value("spherical_position_type", tcc::PositionElementTypes::spherical_position)
        .value("geodetic_position_type", tcc::PositionElementTypes::geodetic_position);

    // ========================================================================
    // Main Conversion Functions
    // ========================================================================

    // Keplerian <-> Cartesian
    function("astro_element_conversion_keplerian_to_cartesian", &keplerianToCartesian);
    function("astro_element_conversion_cartesian_to_keplerian", &cartesianToKeplerian);

    // Spherical <-> Cartesian (state - 6 elements)
    function("astro_element_conversion_spherical_to_cartesian_state", &sphericalToCartesianState);
    function("astro_element_conversion_cartesian_to_spherical_state", &cartesianToSphericalState);

    // MEE conversions
    function("astro_element_conversion_mee_to_cartesian", &meeToCartesian);
    function("astro_element_conversion_cartesian_to_mee", &cartesianToMee);
    function("astro_element_conversion_keplerian_to_mee", &keplerianToMee);
    function("astro_element_conversion_mee_to_keplerian", &meeToKeplerian);

    // Geodetic coordinate conversions
    function("astro_element_conversion_cartesian_to_geodetic", &cartesianToGeodetic);
    function("astro_element_conversion_geodetic_to_cartesian", &geodeticToCartesian);

    // Spherical position conversions (3 elements)
    function("astro_element_conversion_cartesian_to_spherical_position", &cartesianToSphericalPosition);
    function("astro_element_conversion_spherical_to_cartesian_position", &sphericalToCartesianPosition);

    // ========================================================================
    // Anomaly Conversion Functions
    // ========================================================================

    // Elliptical orbit anomaly conversions
    function("astro_element_conversion_mean_to_true_anomaly", &meanToTrueAnomaly);
    function("astro_element_conversion_true_to_mean_anomaly", &trueToMeanAnomaly);
    function("astro_element_conversion_mean_to_eccentric_anomaly", &meanToEccentricAnomaly);
    function("astro_element_conversion_eccentric_to_mean_anomaly", &eccentricToMeanAnomaly);
    function("astro_element_conversion_eccentric_to_true_anomaly", &eccentricToTrueAnomaly);
    function("astro_element_conversion_true_to_eccentric_anomaly", &trueToEccentricAnomaly);

    // Hyperbolic orbit anomaly conversions
    function("astro_element_conversion_mean_to_hyperbolic_eccentric_anomaly", &meanToHyperbolicEccentricAnomaly);
    function("astro_element_conversion_hyperbolic_eccentric_to_mean_anomaly", &hyperbolicEccentricToMeanAnomaly);
    function("astro_element_conversion_hyperbolic_eccentric_to_true_anomaly", &hyperbolicEccentricToTrueAnomaly);
    function("astro_element_conversion_true_to_hyperbolic_eccentric_anomaly", &trueToHyperbolicEccentricAnomaly);

    // ========================================================================
    // USM Element Conversions
    // ========================================================================
    function("astro_element_conversion_cartesian_to_usm_em", &cartesianToUsmEm);
    function("astro_element_conversion_cartesian_to_usm_7", &cartesianToUsm7);
    function("astro_element_conversion_cartesian_to_usm_6", &cartesianToUsm6);
    function("astro_element_conversion_usm_em_to_cartesian", &usmEmToCartesian);
    function("astro_element_conversion_usm_7_to_cartesian", &usm7ToCartesian);
    function("astro_element_conversion_usm_6_to_cartesian", &usm6ToCartesian);

    // ========================================================================
    // Mean Motion and Elapsed Time Conversions
    // ========================================================================
    function("astro_element_conversion_elapsed_time_to_delta_mean_anomaly", &elapsedTimeToDeltaMeanAnomaly);
    function("astro_element_conversion_delta_mean_anomaly_to_elapsed_time", &deltaMeanAnomalyToElapsedTime);
    function("astro_element_conversion_mean_motion_to_semi_major_axis", &meanMotionToSemiMajorAxis);
    function("astro_element_conversion_semi_major_axis_to_mean_motion", &semiMajorAxisToMeanMotion);

    // ========================================================================
    // Quaternion and Rotation Conversions
    // ========================================================================
    function("astro_element_conversion_quaternion_entries_to_rotation_matrix", &quaternionEntriesToRotationMatrix);
    function("astro_element_conversion_rotation_matrix_to_quaternion_entries", &rotationMatrixToQuaternionEntries);
    function("astro_element_conversion_quaternion_to_modified_rodrigues_parameters", &quaternionToModifiedRodriguesParameters);
    function("astro_element_conversion_modified_rodrigues_parameters_to_quaternion", &modifiedRodriguesParametersToQuaternion);
    function("astro_element_conversion_quaternion_to_exponential_map", &quaternionToExponentialMap);
    function("astro_element_conversion_exponential_map_to_quaternion", &exponentialMapToQuaternion);

    // ========================================================================
    // Frame Rotation Conversions
    // ========================================================================
    function("astro_element_conversion_teme_to_j2000", &temeToJ2000);
    function("astro_element_conversion_j2000_to_teme", &j2000ToTeme);
    function("astro_element_conversion_j2000_to_eclipj2000", &j2000ToEclipJ2000);
    function("astro_element_conversion_eclipj2000_to_j2000", &eclipJ2000ToJ2000);

    // ========================================================================
    // Elementwise Conversion Functions
    // ========================================================================
    function("astro_element_conversion_keplerian_to_cartesian_elementwise", &keplerianToCartesianElementwise);
    function("astro_element_conversion_spherical_to_cartesian_elementwise", &sphericalToCartesianElementwise);

    // ========================================================================
    // Spherical Orbital Conversions
    // ========================================================================
    function("astro_element_conversion_spherical_orbital_to_cartesian", &sphericalOrbitalToCartesian);
    function("astro_element_conversion_cartesian_to_spherical_orbital", &cartesianToSphericalOrbital);

    // ========================================================================
    // MEE with Auto Singularity Detection
    // ========================================================================
    function("astro_element_conversion_flip_mee_singularity", &flipMeeSingularity);
    function("astro_element_conversion_keplerian_to_mee_auto", &keplerianToMeeAuto);
    function("astro_element_conversion_cartesian_to_mee_auto", &cartesianToMeeAuto);
}

#endif // __EMSCRIPTEN__
