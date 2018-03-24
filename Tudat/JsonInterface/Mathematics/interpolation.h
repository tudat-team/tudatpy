/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include "Tudat/SimulationSetup/EnvironmentSetup/createGravityFieldVariations.h"
#include "Tudat/JsonInterface/Support/valueAccess.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace interpolators
{

// OneDimensionalInterpolatorTypes

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
inline void to_json( nlohmann::json& jsonObject, const OneDimensionalInterpolatorTypes& oneDimensionalInterpolatorType )
{
    jsonObject = json_interface::stringFromEnum( oneDimensionalInterpolatorType, oneDimensionalInterpolatorTypes );
}

//! Convert `json` to `OneDimensionalInterpolatorTypes`.
inline void from_json( const nlohmann::json& jsonObject, OneDimensionalInterpolatorTypes& oneDimensionalInterpolatorType )
{
    oneDimensionalInterpolatorType =
            json_interface::enumFromString( jsonObject, oneDimensionalInterpolatorTypes );
}


// AvailableLookupScheme

//! Map of `AvailableLookupScheme`s string representations.
static std::map< AvailableLookupScheme, std::string > lookupSchemeTypes =
{
    { huntingAlgorithm, "huntingAlgorithm" },
    { binarySearch, "binarySearch" }
};

//! `AvailableLookupScheme`s not supported by `json_interface`.
static std::vector< AvailableLookupScheme > unsupportedLookupSchemeTypes = { };

//! Convert `AvailableLookupScheme` to `json`.
inline void to_json( nlohmann::json& jsonObject, const AvailableLookupScheme& availableLookupScheme )
{
    jsonObject = json_interface::stringFromEnum( availableLookupScheme, lookupSchemeTypes );
}

//! Convert `json` to `AvailableLookupScheme`.
inline void from_json( const nlohmann::json& jsonObject, AvailableLookupScheme& availableLookupScheme )
{
    availableLookupScheme = json_interface::enumFromString( jsonObject, lookupSchemeTypes );
}


// LagrangeInterpolatorBoundaryHandling

//! Map of `LagrangeInterpolatorBoundaryHandling`s string representations.
static std::map< LagrangeInterpolatorBoundaryHandling, std::string > lagrangeInterpolatorBoundaryHandlings =
{
    { lagrange_cubic_spline_boundary_interpolation, "cubicSplineBoundary" },
    { lagrange_no_boundary_interpolation, "noBoundary" }
};

//! `LagrangeInterpolatorBoundaryHandling`s not supported by `json_interface`.
static std::vector< LagrangeInterpolatorBoundaryHandling > unsupportedLagrangeInterpolatorBoundaryHandlings = { };

//! Convert `LagrangeInterpolatorBoundaryHandling` to `json`.
inline void to_json( nlohmann::json& jsonObject,
                     const LagrangeInterpolatorBoundaryHandling& lagrangeInterpolatorBoundaryHandling )
{
    jsonObject = json_interface::stringFromEnum( lagrangeInterpolatorBoundaryHandling,
                                                 lagrangeInterpolatorBoundaryHandlings );
}

//! Convert `json` to `LagrangeInterpolatorBoundaryHandling`.
inline void from_json( const nlohmann::json& jsonObject,
                       LagrangeInterpolatorBoundaryHandling& lagrangeInterpolatorBoundaryHandling )
{
    lagrangeInterpolatorBoundaryHandling =
            json_interface::enumFromString( jsonObject, lagrangeInterpolatorBoundaryHandlings );
}


// DataMapSettings

//! Create a `json` object from a shared pointer to a `DataMapSettings` object.
template< typename I, typename D >
void to_json( nlohmann::json& jsonObject, const boost::shared_ptr< DataMapSettings< I, D > >& dataMapSettings )
{
    if ( ! dataMapSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Interpolation::DataMap;

    boost::shared_ptr< HermiteDataSettings< I, D > > hermiteDataSettings =
            boost::dynamic_pointer_cast< HermiteDataSettings< I, D > >( dataMapSettings );
    if ( hermiteDataSettings )
    {
        jsonObject[ K::map ] = hermiteDataSettings->getDataMap( );
        jsonObject[ K::dependentVariableFirstDerivativeValues ] =
                hermiteDataSettings->firstDerivativeOfDependentVariables_;
        return;
    }

    boost::shared_ptr< FromFileDataMapSettings< D > > fromFileDataMapSettings =
            boost::dynamic_pointer_cast< FromFileDataMapSettings< D > >( dataMapSettings );
    if ( fromFileDataMapSettings )
    {
        jsonObject[ K::file ] = boost::filesystem::path( fromFileDataMapSettings->relativeFilePath_ );
        return;
    }

    boost::shared_ptr< IndependentDependentDataMapSettings< I, D > > independentDependentDataMapSettings =
            boost::dynamic_pointer_cast< IndependentDependentDataMapSettings< I, D > >( dataMapSettings );
    if ( independentDependentDataMapSettings )
    {
        jsonObject[ K::independentVariableValues ] = independentDependentDataMapSettings->independentVariableValues_;
        jsonObject[ K::dependentVariableValues ] = independentDependentDataMapSettings->dependentVariableValues_;
        return;
    }

    jsonObject[ K::map ] = dataMapSettings->getDataMap( );
}

//! Create a shared pointer to a `DataMapSettings` object from a `json` object.
template< typename I, typename D >
void from_json( const nlohmann::json& jsonObject, boost::shared_ptr< DataMapSettings< I, D > >& dataMapSettings )
{
    using namespace json_interface;
    using K = Keys::Interpolation::DataMap;

    if ( isDefined( jsonObject, K::dependentVariableFirstDerivativeValues ) )
    {
        dataMapSettings = boost::make_shared< HermiteDataSettings< I, D > >(
                    getValue< std::map< I, D > >( jsonObject, K::map ),
                    getValue< std::vector< D > >( jsonObject, K::dependentVariableFirstDerivativeValues ) );
        return;
    }

    if ( isDefined( jsonObject, K::file ) )
    {
        dataMapSettings = boost::make_shared< FromFileDataMapSettings< D > >(
                    getValue< boost::filesystem::path >( jsonObject, K::file ).string( ) );
        return;
    }

    if ( isDefined( jsonObject, K::independentVariableValues ) || isDefined( jsonObject, K::dependentVariableValues ) )
    {
        dataMapSettings = boost::make_shared< IndependentDependentDataMapSettings< I, D > >(
                    getValue< std::vector< I > >( jsonObject, K::independentVariableValues ),
                    getValue< std::vector< D > >( jsonObject, K::dependentVariableValues ) );
        return;
    }

    dataMapSettings = boost::make_shared< DataMapSettings< I, D > >(
                getValue< std::map< I, D > >( jsonObject, K::map ) );
    return;
}


// InterpolatorSettings

//! Create a `json` object from a shared pointer to a `InterpolatorSettings` object.
void to_json( nlohmann::json& jsonObject, const boost::shared_ptr< InterpolatorSettings >& interpolatorSettings );

//! Create a shared pointer to a `InterpolatorSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, boost::shared_ptr< InterpolatorSettings >& interpolatorSettings );


// DataInterpolationSettings

//! Create a `json` object from a shared pointer to a `DataInterpolationSettings` object.
template< typename I, typename D >
void to_json( nlohmann::json& jsonObject,
              const boost::shared_ptr< DataInterpolationSettings< I, D > >& dataInterpolationSettings )
{
    if ( ! dataInterpolationSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Interpolation::DataInterpolation;

    jsonObject[ K::data ] = dataInterpolationSettings->dataSettings_;
    jsonObject[ K::interpolator ] = dataInterpolationSettings->interpolatorSettings_;
}

//! Create a shared pointer to a `DataInterpolationSettings` object from a `json` object.
template< typename I, typename D >
void from_json( const nlohmann::json& jsonObject,
                boost::shared_ptr< DataInterpolationSettings< I, D > >& dataInterpolationSettings )
{
    using namespace json_interface;
    using K = Keys::Interpolation::DataInterpolation;

    dataInterpolationSettings = boost::make_shared< DataInterpolationSettings< I, D > >(
                getValue< boost::shared_ptr< DataMapSettings< I, D > > >( jsonObject, K::data ),
                getValue< boost::shared_ptr< InterpolatorSettings > >( jsonObject, K::interpolator ) );
}

} // namespace interpolators


namespace simulation_setup
{

// ModelInterpolationSettings

//! Create a `json` object from a shared pointer to a `ModelInterpolationSettings` object.
void to_json( nlohmann::json& jsonObject, const boost::shared_ptr< ModelInterpolationSettings >& modelInterpolationSettings );

//! Create a shared pointer to a `ModelInterpolationSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, boost::shared_ptr< ModelInterpolationSettings >& modelInterpolationSettings );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_INTERPOLATOR_H
