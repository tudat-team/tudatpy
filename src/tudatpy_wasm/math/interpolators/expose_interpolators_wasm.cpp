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
#include "../../shared_ptr_wasm.h"

#include <tudat/math/interpolators.h>

namespace ti = tudat::interpolators;

WASM_MODULE_PATH("math_interpolators")

EMSCRIPTEN_BINDINGS(tudatpy_math_interpolators) {
    using namespace emscripten;

    // Enums
    enum_<ti::BoundaryInterpolationType>("math_interpolators_BoundaryInterpolationType")
        .value("throw_exception_at_boundary", ti::throw_exception_at_boundary)
        .value("use_boundary_value", ti::use_boundary_value)
        .value("use_boundary_value_with_warning", ti::use_boundary_value_with_warning)
        .value("extrapolate_at_boundary", ti::extrapolate_at_boundary)
        .value("extrapolate_at_boundary_with_warning", ti::extrapolate_at_boundary_with_warning);

    enum_<ti::AvailableLookupScheme>("math_interpolators_AvailableLookupScheme")
        .value("hunting_algorithm", ti::huntingAlgorithm)
        .value("binary_search", ti::binarySearch);

    enum_<ti::LagrangeInterpolatorBoundaryHandling>("math_interpolators_LagrangeInterpolatorBoundaryHandling")
        .value("lagrange_cubic_spline_boundary_interpolation", ti::lagrange_cubic_spline_boundary_interpolation)
        .value("lagrange_no_boundary_interpolation", ti::lagrange_no_boundary_interpolation);

    // Base InterpolatorSettings class
    class_<ti::InterpolatorSettings>("math_interpolators_InterpolatorSettings")
        .smart_ptr<std::shared_ptr<ti::InterpolatorSettings>>("shared_ptr_InterpolatorSettings");

    // LagrangeInterpolatorSettings
    class_<ti::LagrangeInterpolatorSettings, base<ti::InterpolatorSettings>>("math_interpolators_LagrangeInterpolatorSettings")
        .smart_ptr<std::shared_ptr<ti::LagrangeInterpolatorSettings>>("shared_ptr_LagrangeInterpolatorSettings");

    // Factory functions (non-templated inline functions)
    function("math_interpolators_linear_interpolation",
        &ti::linearInterpolation);

    function("math_interpolators_cubic_spline_interpolation",
        &ti::cubicSplineInterpolation);

    function("math_interpolators_piecewise_constant_interpolation",
        &ti::piecewiseConstantInterpolation);

    function("math_interpolators_lagrange_interpolation",
        &ti::lagrangeInterpolation);

    function("math_interpolators_hermite_interpolation",
        &ti::hermiteInterpolation);
}

#endif
