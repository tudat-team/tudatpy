/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_INTERPOLATOR_H
#define TUDAT_JSONINTERFACE_INTERPOLATOR_H

#include <Tudat/SimulationSetup/EnvironmentSetup/createGravityFieldVariations.h>

#include "Tudat/External/JsonInterface/Support/valueAccess.h"
#include "Tudat/External/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace interpolators
{

//! Map of `OneDimensionalInterpolatorTypes` supported by `json_interface`.
static std::map< std::string, OneDimensionalInterpolatorTypes > oneDimensionalInterpolatorTypes =
{
    { "linear",            linear_interpolator },
    { "cubicSpline",       cubic_spline_interpolator },
    { "lagrange",          lagrange_interpolator },
    { "hermiteSpline",     hermite_spline_interpolator },
    { "piecewiseConstant", piecewise_constant_interpolator }
};

//! Convert `OneDimensionalInterpolatorTypes` to `json`.
void to_json( json& jsonObject, const OneDimensionalInterpolatorTypes& oneDimensionalInterpolatorType );

//! Convert `json` to `OneDimensionalInterpolatorTypes`.
void from_json( const json& jsonObject, OneDimensionalInterpolatorTypes& oneDimensionalInterpolatorType );


//! Map of `AvailableLookupScheme`s supported by `json_interface`.
static std::map< std::string, AvailableLookupScheme > availableLookupSchemes =
{
    { "huntingAlgorithm", huntingAlgorithm },
    { "binarySearch",     binarySearch }
};

//! Convert `AvailableLookupScheme` to `json`.
void to_json( json& jsonObject, const AvailableLookupScheme& availableLookupScheme );

//! Convert `json` to `AvailableLookupScheme`.
void from_json( const json& jsonObject, AvailableLookupScheme& availableLookupScheme );


//! Map of `LagrangeInterpolatorBoundaryHandling`s supported by `json_interface`.
static std::map< std::string, LagrangeInterpolatorBoundaryHandling > lagrangeInterpolatorBoundaryHandlings =
{
    { "cubicSplineBoundary", lagrange_cubic_spline_boundary_interpolation },
    { "noBoundary",          lagrange_no_boundary_interpolation }
};

//! Convert `LagrangeInterpolatorBoundaryHandling` to `json`.
void to_json( json& jsonObject, const LagrangeInterpolatorBoundaryHandling& lagrangeInterpolatorBoundaryHandling );

//! Convert `json` to `LagrangeInterpolatorBoundaryHandling`.
void from_json( const json& jsonObject, LagrangeInterpolatorBoundaryHandling& lagrangeInterpolatorBoundaryHandling );


//! Create a `json` object from a shared pointer to a `InterpolatorSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< InterpolatorSettings >& interpolatorSettings );

//! Create a shared pointer to a `InterpolatorSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< InterpolatorSettings >& interpolatorSettings );

} // namespace interpolators


namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `ModelInterpolationSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< ModelInterpolationSettings >& modelInterpolationSettings );

//! Create a shared pointer to a `ModelInterpolationSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< ModelInterpolationSettings >& modelInterpolationSettings );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_INTERPOLATOR_H
