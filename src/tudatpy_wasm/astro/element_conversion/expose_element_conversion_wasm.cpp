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

namespace toec = tudat::orbital_element_conversions;
namespace tcc = tudat::coordinate_conversions;
namespace tla = tudat::linear_algebra;
namespace te = tudat::ephemerides;
namespace tba = tudat::basic_astrodynamics;

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

// Note: elapsedTimeFromMeanMotionAndAnomaly removed - API mismatch

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
}

#endif // __EMSCRIPTEN__
