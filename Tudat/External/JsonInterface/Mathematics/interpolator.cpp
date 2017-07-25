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
        using Keys = Keys::Interpolator;

        // Common properties
        jsonObject[ Keys::type ] = interpolatorSettings->getInterpolatorType( );
        jsonObject[ Keys::lookupScheme ] = interpolatorSettings->getSelectedLookupScheme( );
        jsonObject[ Keys::useLongDoubleTimeStep ] = interpolatorSettings->getUseLongDoubleTimeStep( );

        /// LagrangeInterpolatorSettings
        boost::shared_ptr< LagrangeInterpolatorSettings > lagrangeInterpolatorSettings =
                boost::dynamic_pointer_cast< LagrangeInterpolatorSettings >( interpolatorSettings );
        if ( lagrangeInterpolatorSettings )
        {
            jsonObject[ Keys::order ] = lagrangeInterpolatorSettings->getInterpolatorOrder( );
            jsonObject[ Keys::boundaryHandling ] = lagrangeInterpolatorSettings->getBoundaryHandling( );
            return;
        }
    }
}

/*
//! Convert `json` to `InterpolatorSettings` shared pointer.
void from_json( const json& jsonObject, boost::shared_ptr< InterpolatorSettings >& interpolatorSettings )
{
    interpolatorSettings = json_interface::createInterpolatorSettings( jsonObject );
}
*/

} // namespace interpolators



namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `ModelInterpolationSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< ModelInterpolationSettings >& modelInterpolationSettings )
{
    if ( modelInterpolationSettings )
    {
        using namespace json_interface;
        using Keys = Keys::ModelInterpolation;

        jsonObject[ Keys::initialTime ] = modelInterpolationSettings->initialTime_;
        jsonObject[ Keys::finalTime ] = modelInterpolationSettings->finalTime_;
        jsonObject[ Keys::timeStep ] = modelInterpolationSettings->timeStep_;
        jsonObject[ Keys::interpolator ] = modelInterpolationSettings->interpolatorSettings_;
    }
}

/*
//! Convert `json` to `ModelInterpolationSettings` shared pointer.
void from_json( const json& jsonObject, boost::shared_ptr< ModelInterpolationSettings >& modelInterpolationSettings )
{
    modelInterpolationSettings = json_interface::createModelInterpolationSettings( jsonObject );
}
*/

} // namespace simulation_setup



namespace json_interface
{

//! Create a shared pointer to a `InterpolatorSettings` object from a `json` object.
boost::shared_ptr< interpolators::InterpolatorSettings > createInterpolatorSettings(
        const json& settings, const KeyTree& keyTree,
        const boost::shared_ptr< interpolators::InterpolatorSettings >& fallback )
{
    if ( ! defined( settings, keyTree ) )
    {
        return fallback;
    }
    else
    {
        using namespace interpolators;
        using Keys = Keys::Interpolator;

        // Get interpolator type
        const OneDimensionalInterpolatorTypes oneDimensionalInterpolatorType =
                getValue< OneDimensionalInterpolatorTypes >( settings, keyTree + Keys::type );

        switch ( oneDimensionalInterpolatorType ) {
        case linear_interpolator:
        case cubic_spline_interpolator:
        case hermite_spline_interpolator:
        case piecewise_constant_interpolator:
        {
            InterpolatorSettings defaults( linear_interpolator );
            return boost::make_shared< InterpolatorSettings >(
                        oneDimensionalInterpolatorType,
                        getValue( settings, keyTree + Keys::lookupScheme,
                                  defaults.getSelectedLookupScheme( ) ),
                        getValue( settings, keyTree + Keys::useLongDoubleTimeStep,
                                  defaults.getUseLongDoubleTimeStep( ) ) );
        }
        case lagrange_interpolator:
        {
            LagrangeInterpolatorSettings defaults( 0 );
            return boost::make_shared< LagrangeInterpolatorSettings >(
                        getValue< double >( settings, keyTree + Keys::order ),
                        getValue( settings, keyTree + Keys::useLongDoubleTimeStep,
                                  defaults.getUseLongDoubleTimeStep( ) ),
                        getValue( settings, keyTree + Keys::lookupScheme,
                                  defaults.getSelectedLookupScheme( ) ),
                        getValue( settings, keyTree + Keys::boundaryHandling,
                                  defaults.getBoundaryHandling( ) ) );
        }
        default:
            throw std::runtime_error( stringFromEnum( oneDimensionalInterpolatorType, oneDimensionalInterpolatorTypes )
                                      + " not supported by json_interface." );
        }
    }
}

//! Create a shared pointer to a `ModelInterpolationSettings` object from a `json` object.
boost::shared_ptr< simulation_setup::ModelInterpolationSettings > createModelInterpolationSettings(
        const json& settings, const KeyTree& keyTree,
        const boost::shared_ptr< simulation_setup::ModelInterpolationSettings >& fallback )
{
    if ( ! defined( settings, keyTree ) )
    {
        return fallback;
    }
    else
    {
        using namespace simulation_setup;
        using Keys = Keys::ModelInterpolation;

        ModelInterpolationSettings defaults;
        return boost::make_shared< ModelInterpolationSettings >(
                    getEpoch( settings, keyTree + Keys::initialTime, defaults.initialTime_ ),
                    getEpoch( settings, keyTree + Keys::finalTime, defaults.finalTime_ ),
                    getNumeric( settings, keyTree + Keys::timeStep, defaults.timeStep_ ),
                    createInterpolatorSettings( settings, keyTree + Keys::interpolator,
                                                defaults.interpolatorSettings_ ) );
    }
}

} // namespace json_interface

} // namespace tudat
