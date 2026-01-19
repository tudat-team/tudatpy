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

#include <tudat/math/basic/legendrePolynomials.h>
#include <tudat/astro/gravitation/sphericalHarmonicsGravityField.h>

namespace tbm = tudat::basic_mathematics;
namespace tg = tudat::gravitation;

WASM_MODULE_PATH("astro_gravitation")

namespace {
using namespace tudatpy_wasm;

// Wrapper for normalize_spherical_harmonic_coefficients
// Returns a pair of MatrixXdWrapper (cosine, sine coefficients)
emscripten::val normalizeSphericalHarmonicCoefficients(
    const MatrixXdWrapper& unnormalizedCosine,
    const MatrixXdWrapper& unnormalizedSine) {

    auto result = tbm::convertUnnormalizedToGeodesyNormalizedCoefficients(
        unnormalizedCosine.data, unnormalizedSine.data);

    emscripten::val output = emscripten::val::object();
    output.set("cosine", MatrixXdWrapper(std::get<0>(result)));
    output.set("sine", MatrixXdWrapper(std::get<1>(result)));
    return output;
}

// Wrapper for unnormalize_spherical_harmonic_coefficients
emscripten::val unnormalizeSphericalHarmonicCoefficients(
    const MatrixXdWrapper& normalizedCosine,
    const MatrixXdWrapper& normalizedSine) {

    auto result = tbm::convertGeodesyNormalizedToUnnormalizedCoefficients(
        normalizedCosine.data, normalizedSine.data);

    emscripten::val output = emscripten::val::object();
    output.set("cosine", MatrixXdWrapper(std::get<0>(result)));
    output.set("sine", MatrixXdWrapper(std::get<1>(result)));
    return output;
}

// Wrapper for spherical_harmonic_coefficients_from_inertia
emscripten::val sphericalHarmonicCoefficientsFromInertia(
    const Matrix3dWrapper& inertiaTensor,
    double gravitationalParameter,
    double referenceRadius,
    bool outputNormalizedCoefficients) {

    auto result = tg::getDegreeTwoSphericalHarmonicCoefficients(
        inertiaTensor.data,
        gravitationalParameter,
        referenceRadius,
        2,
        outputNormalizedCoefficients);

    emscripten::val output = emscripten::val::object();
    output.set("cosine", MatrixXdWrapper(std::get<0>(result)));
    output.set("sine", MatrixXdWrapper(std::get<1>(result)));
    output.set("meanMomentOfInertia", std::get<2>(result));
    return output;
}

}

EMSCRIPTEN_BINDINGS(tudatpy_astro_gravitation) {
    using namespace emscripten;

    function("astro_gravitation_legendre_normalization_factor",
        &tbm::calculateLegendreGeodesyNormalizationFactor);

    function("astro_gravitation_normalize_spherical_harmonic_coefficients",
        &normalizeSphericalHarmonicCoefficients);

    function("astro_gravitation_unnormalize_spherical_harmonic_coefficients",
        &unnormalizeSphericalHarmonicCoefficients);

    function("astro_gravitation_spherical_harmonic_coefficients_from_inertia",
        &sphericalHarmonicCoefficientsFromInertia);
}

#endif
