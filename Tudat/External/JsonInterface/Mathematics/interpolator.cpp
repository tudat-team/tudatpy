/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "interpolator.h"

namespace tudat
{

namespace interpolators
{

//! Convert `OneDimensionalInterpolatorTypes` to `json`.
void to_json( json& jsonObject, const OneDimensionalInterpolatorTypes& oneDimensionalInterpolatorType )
{
    jsonObject = json_interface::stringFromEnum( oneDimensionalInterpolatorType, oneDimensionalInterpolatorTypes );
}

//! Convert `json` to `OneDimensionalInterpolatorTypes`.
void from_json( const json& jsonObject, OneDimensionalInterpolatorTypes& oneDimensionalInterpolatorType )
{
    oneDimensionalInterpolatorType =
            json_interface::enumFromString( jsonObject.get< std::string >( ), oneDimensionalInterpolatorTypes );
}


//! Convert `AvailableLookupScheme` to `json`.
void to_json( json& jsonObject, const AvailableLookupScheme& availableLookupScheme )
{
    jsonObject = json_interface::stringFromEnum( availableLookupScheme, availableLookupSchemes );
}

//! Convert `json` to `AvailableLookupScheme`.
void from_json( const json& jsonObject, AvailableLookupScheme& availableLookupScheme )
{
    availableLookupScheme = json_interface::enumFromString( jsonObject.get< std::string >( ), availableLookupSchemes );
}


//! Convert `LagrangeInterpolatorBoundaryHandling` to `json`.
void to_json( json& jsonObject, const LagrangeInterpolatorBoundaryHandling& lagrangeInterpolatorBoundaryHandling )
{
    jsonObject = json_interface::stringFromEnum( lagrangeInterpolatorBoundaryHandling,
                                                 lagrangeInterpolatorBoundaryHandlings );
}

//! Convert `json` to `AvailableLookupScheme`.
void from_json( const json& jsonObject, LagrangeInterpolatorBoundaryHandling& lagrangeInterpolatorBoundaryHandling )
{
    lagrangeInterpolatorBoundaryHandling =
            json_interface::enumFromString( jsonObject.get< std::string >( ), lagrangeInterpolatorBoundaryHandlings );
}


//! Create a `json` object from a shared pointer to a `InterpolatorSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< InterpolatorSettings >& interpolatorSettings )
{
    if ( interpolatorSettings )
    {
        using namespace json_interface;
        using K = Keys::Interpolator;

        // Common properties
        jsonObject[ K::type ] = interpolatorSettings->getInterpolatorType( );
        jsonObject[ K::lookupScheme ] = interpolatorSettings->getSelectedLookupScheme( );
        jsonObject[ K::useLongDoubleTimeStep ] = interpolatorSettings->getUseLongDoubleTimeStep( );

        /// LagrangeInterpolatorSettings
        boost::shared_ptr< LagrangeInterpolatorSettings > lagrangeInterpolatorSettings =
                boost::dynamic_pointer_cast< LagrangeInterpolatorSettings >( interpolatorSettings );
        if ( lagrangeInterpolatorSettings )
        {
            jsonObject[ K::order ] = lagrangeInterpolatorSettings->getInterpolatorOrder( );
            jsonObject[ K::boundaryHandling ] = lagrangeInterpolatorSettings->getBoundaryHandling( );
            return;
        }
    }
}

//! Create a shared pointer to a `InterpolatorSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< InterpolatorSettings >& interpolatorSettings )
{
    using namespace json_interface;
    using K = Keys::Interpolator;

    // Get interpolator type
    const OneDimensionalInterpolatorTypes oneDimensionalInterpolatorType =
            getValue< OneDimensionalInterpolatorTypes >( jsonObject, K::type );

    switch ( oneDimensionalInterpolatorType ) {
    case linear_interpolator:
    case cubic_spline_interpolator:
    case hermite_spline_interpolator:
    case piecewise_constant_interpolator:
    {
        InterpolatorSettings defaults( linear_interpolator );
        interpolatorSettings = boost::make_shared< InterpolatorSettings >(
                    oneDimensionalInterpolatorType,
                    getValue( jsonObject, K::lookupScheme, defaults.getSelectedLookupScheme( ) ),
                    getValue( jsonObject, K::useLongDoubleTimeStep, defaults.getUseLongDoubleTimeStep( ) ) );
        return;
    }
    case lagrange_interpolator:
    {
        LagrangeInterpolatorSettings defaults( 0 );
        interpolatorSettings = boost::make_shared< LagrangeInterpolatorSettings >(
                    getValue< double >( jsonObject, K::order ),
                    getValue( jsonObject, K::useLongDoubleTimeStep, defaults.getUseLongDoubleTimeStep( ) ),
                    getValue( jsonObject, K::lookupScheme, defaults.getSelectedLookupScheme( ) ),
                    getValue( jsonObject, K::boundaryHandling, defaults.getBoundaryHandling( ) ) );
        return;
    }
    default:
        throw std::runtime_error( stringFromEnum( oneDimensionalInterpolatorType, oneDimensionalInterpolatorTypes )
                                  + " not supported by json_interface." );
    }
}

} // namespace interpolators


namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `ModelInterpolationSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< ModelInterpolationSettings >& modelInterpolationSettings )
{
    if ( modelInterpolationSettings )
    {
        using namespace json_interface;
        using K = Keys::ModelInterpolation;

        jsonObject[ K::initialTime ] = modelInterpolationSettings->initialTime_;
        jsonObject[ K::finalTime ] = modelInterpolationSettings->finalTime_;
        jsonObject[ K::timeStep ] = modelInterpolationSettings->timeStep_;
        jsonObject[ K::interpolator ] = modelInterpolationSettings->interpolatorSettings_;
    }
}

//! Create a shared pointer to a `ModelInterpolationSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< ModelInterpolationSettings >& modelInterpolationSettings )
{
    using namespace json_interface;
    using K = Keys::ModelInterpolation;

    ModelInterpolationSettings defaults;
    modelInterpolationSettings = boost::make_shared< ModelInterpolationSettings >(
                getEpoch( jsonObject, K::initialTime, defaults.initialTime_ ),
                getEpoch( jsonObject, K::finalTime, defaults.finalTime_ ),
                getNumeric( jsonObject, K::timeStep, defaults.timeStep_ ),
                getValue( jsonObject, K::interpolator, defaults.interpolatorSettings_ ) );
    return;
}

} // namespace simulation_setup

} // namespace tudat
