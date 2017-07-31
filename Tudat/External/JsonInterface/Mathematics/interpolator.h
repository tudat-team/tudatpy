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

//! Map of `OneDimensionalInterpolatorTypes` string representations.
static std::map< OneDimensionalInterpolatorTypes, std::string > oneDimensionalInterpolatorTypes =
{
    { linear_interpolator, "linear" },
    { cubic_spline_interpolator, "cubicSpline" },
    { lagrange_interpolator, "lagrange" },
    { hermite_spline_interpolator, "hermiteSpline" },
    { piecewise_constant_interpolator, "piecewiseConstant" }
};

//! `OneDimensionalInterpolatorTypes` not supported by `json_interface`.
static std::vector< OneDimensionalInterpolatorTypes > unsupportedOneDimensionalInterpolatorTypes = { };

//! Convert `OneDimensionalInterpolatorTypes` to `json`.
inline void to_json( json& jsonObject, const OneDimensionalInterpolatorTypes& oneDimensionalInterpolatorType )
{
    jsonObject = json_interface::stringFromEnum( oneDimensionalInterpolatorType, oneDimensionalInterpolatorTypes );
}

//! Convert `json` to `OneDimensionalInterpolatorTypes`.
inline void from_json( const json& jsonObject, OneDimensionalInterpolatorTypes& oneDimensionalInterpolatorType )
{
    oneDimensionalInterpolatorType =
            json_interface::enumFromString( jsonObject.get< std::string >( ), oneDimensionalInterpolatorTypes );
}


//! Map of `AvailableLookupScheme`s string representations.
static std::map< AvailableLookupScheme, std::string > lookupSchemeTypes =
{
    { huntingAlgorithm, "huntingAlgorithm" },
    { binarySearch, "binarySearch" }
};

//! `AvailableLookupScheme`s not supported by `json_interface`.
static std::vector< AvailableLookupScheme > unsupportedLookupSchemeTypes = { };

//! Convert `AvailableLookupScheme` to `json`.
inline void to_json( json& jsonObject, const AvailableLookupScheme& availableLookupScheme )
{
    jsonObject = json_interface::stringFromEnum( availableLookupScheme, lookupSchemeTypes );
}

//! Convert `json` to `AvailableLookupScheme`.
inline void from_json( const json& jsonObject, AvailableLookupScheme& availableLookupScheme )
{
    availableLookupScheme = json_interface::enumFromString( jsonObject.get< std::string >( ), lookupSchemeTypes );
}


//! Map of `LagrangeInterpolatorBoundaryHandling`s string representations.
static std::map< LagrangeInterpolatorBoundaryHandling, std::string > lagrangeInterpolatorBoundaryHandlings =
{
    { lagrange_cubic_spline_boundary_interpolation, "cubicSplineBoundary" },
    { lagrange_no_boundary_interpolation, "noBoundary" }
};

//! Convert `LagrangeInterpolatorBoundaryHandling` to `json`.
inline void to_json( json& jsonObject,
                     const LagrangeInterpolatorBoundaryHandling& lagrangeInterpolatorBoundaryHandling )
{
    jsonObject = json_interface::stringFromEnum( lagrangeInterpolatorBoundaryHandling,
                                                 lagrangeInterpolatorBoundaryHandlings );
}

//! Convert `json` to `LagrangeInterpolatorBoundaryHandling`.
inline void from_json( const json& jsonObject,
                       LagrangeInterpolatorBoundaryHandling& lagrangeInterpolatorBoundaryHandling )
{
    lagrangeInterpolatorBoundaryHandling =
            json_interface::enumFromString( jsonObject.get< std::string >( ), lagrangeInterpolatorBoundaryHandlings );
}


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
