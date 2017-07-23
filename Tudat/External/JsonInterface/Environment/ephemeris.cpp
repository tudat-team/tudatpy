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

#include "ephemeris.h"

#include "Tudat/External/JsonInterface/Mathematics/interpolator.h"

namespace tudat
{

namespace ephemerides
{

//! Convert `ApproximatePlanetPositionsBase::BodiesWithEphemerisData` to `json`.
void to_json( json& jsonObject, const ApproximatePlanetPositionsBase::BodiesWithEphemerisData& bodyWithEphemerisData )
{
    jsonObject = json_interface::stringFromEnum( bodyWithEphemerisData, bodiesWithEphemerisData );
}

//! Convert `json` to `ApproximatePlanetPositionsBase::BodiesWithEphemerisData`.
void from_json( const json& jsonObject, ApproximatePlanetPositionsBase::BodiesWithEphemerisData& bodyWithEphemerisData )
{
    bodyWithEphemerisData = json_interface::enumFromString( jsonObject.get< std::string >( ), bodiesWithEphemerisData );
}

} // namespace ephemerides


namespace simulation_setup
{

//! Convert `EphemerisType` to `json`.
void to_json( json& jsonObject, const EphemerisType& ephemerisType )
{
    jsonObject = json( json_interface::stringFromEnum( ephemerisType, ephemerisTypes ) );
}

//! Convert `json` to `EphemerisType`.
void from_json( const json& jsonObject, EphemerisType& ephemerisType )
{
    ephemerisType = json_interface::enumFromString( jsonObject.get< std::string >( ), ephemerisTypes );
}

//! Create a `json` object from a shared pointer to a `EphemerisSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< EphemerisSettings >& ephemerisSettings )
{
    if ( ephemerisSettings )
    {
        using namespace json_interface;
        using Keys = Keys::Body::Ephemeris;

        // Common keys
        jsonObject[ Keys::type ] = ephemerisSettings->getEphemerisType( );
        jsonObject[ Keys::frameOrigin ] = ephemerisSettings->getFrameOrigin( );
        jsonObject[ Keys::frameOrientation ] = ephemerisSettings->getFrameOrientation( );
        jsonObject[ Keys::makeMultiArc ] = ephemerisSettings->getMakeMultiArcEphemeris( );

        /// ApproximatePlanetPositionSettings
        boost::shared_ptr< ApproximatePlanetPositionSettings > approximatePlanetPositionSettings =
                boost::dynamic_pointer_cast< ApproximatePlanetPositionSettings >( ephemerisSettings );
        if ( approximatePlanetPositionSettings )
        {
            jsonObject[ Keys::bodyIdentifier ] = approximatePlanetPositionSettings->getBodyIdentifier( );
            jsonObject[ Keys::useCircularCoplanarApproximation ] =
                    approximatePlanetPositionSettings->getUseCircularCoplanarApproximation( );
            return;
        }

        /// DirectSpiceEphemerisSettings
        boost::shared_ptr< DirectSpiceEphemerisSettings > directSpiceEphemerisSettings =
                boost::dynamic_pointer_cast< DirectSpiceEphemerisSettings >( ephemerisSettings );
        if ( directSpiceEphemerisSettings )
        {
            jsonObject[ Keys::correctForStellarAbberation ] =
                    directSpiceEphemerisSettings->getCorrectForStellarAbberation( );
            jsonObject[ Keys::correctForLightTimeAbberation ] =
                    directSpiceEphemerisSettings->getCorrectForLightTimeAbberation( );
            jsonObject[ Keys::convergeLighTimeAbberation ] =
                    directSpiceEphemerisSettings->getConvergeLighTimeAbberation( );

            /// InterpolatedSpiceEphemerisSettings
            boost::shared_ptr< InterpolatedSpiceEphemerisSettings > interpolatedSpiceEphemerisSettings =
                    boost::dynamic_pointer_cast< InterpolatedSpiceEphemerisSettings >( directSpiceEphemerisSettings );
            if ( interpolatedSpiceEphemerisSettings )
            {
                jsonObject[ Keys::initialTime ] = interpolatedSpiceEphemerisSettings->getInitialTime( );
                jsonObject[ Keys::finalTime ] = interpolatedSpiceEphemerisSettings->getFinalTime( );
                jsonObject[ Keys::timeStep ] = interpolatedSpiceEphemerisSettings->getTimeStep( );
                jsonObject[ Keys::interpolator ] = interpolatedSpiceEphemerisSettings->getInterpolatorSettings( );
                jsonObject[ Keys::useLongDoubleStates ] = interpolatedSpiceEphemerisSettings->getUseLongDoubleStates( );
                return;
            }

            return;
        }

        /// TabulatedEphemerisSettings
        boost::shared_ptr< TabulatedEphemerisSettings > tabulatedEphemerisSettings =
                boost::dynamic_pointer_cast< TabulatedEphemerisSettings >( ephemerisSettings );
        if ( tabulatedEphemerisSettings )
        {
            jsonObject[ Keys::bodyStateHistory ] = tabulatedEphemerisSettings->getBodyStateHistory( );
            jsonObject[ Keys::useLongDoubleStates ] = tabulatedEphemerisSettings->getUseLongDoubleStates( );
            return;
        }

        /// ConstantEphemerisSettings
        boost::shared_ptr< ConstantEphemerisSettings > constantEphemerisSettings =
                boost::dynamic_pointer_cast< ConstantEphemerisSettings >( ephemerisSettings );
        if ( constantEphemerisSettings )
        {
            jsonObject[ Keys::constantState ] = constantEphemerisSettings->getConstantState( );
            return;
        }

        /// KeplerEphemerisSettings
        boost::shared_ptr< KeplerEphemerisSettings > keplerEphemerisSettings =
                boost::dynamic_pointer_cast< KeplerEphemerisSettings >( ephemerisSettings );
        if ( keplerEphemerisSettings )
        {
            jsonObject[ Keys::initialStateInKeplerianElements ] =
                    keplerEphemerisSettings->getInitialStateInKeplerianElements( );
            jsonObject[ Keys::epochOfInitialState ] =
                    keplerEphemerisSettings->getEpochOfInitialState( );
            jsonObject[ Keys::centralBodyGravitationalParameter ] =
                    keplerEphemerisSettings->getCentralBodyGravitationalParameter( );
            jsonObject[ Keys::rootFinderAbsoluteTolerance ] =
                    keplerEphemerisSettings->getRootFinderAbsoluteTolerance( );
            jsonObject[ Keys::rootFinderMaximumNumberOfIterations ] =
                    keplerEphemerisSettings->getRootFinderMaximumNumberOfIterations( );
            return;
        }

        /// CustomEphemerisSettings
        boost::shared_ptr< CustomEphemerisSettings > customEphemerisSettings =
                boost::dynamic_pointer_cast< CustomEphemerisSettings >( ephemerisSettings );
        if ( customEphemerisSettings )
        {
            throw std::runtime_error( "CustomEphemerisSettings not supported by json_interface." ); // FIXME
        }
    }
}

} // namespace simulation_setup


namespace json_interface
{

//! Create a shared pointer to a `EphemerisSettings` object from a `json` object.
boost::shared_ptr< simulation_setup::EphemerisSettings > createEphemerisSettings(
        const json& settings, const KeyTree& keyTree )
{
    using namespace simulation_setup;
    using Keys = Keys::Body::Ephemeris;

    // Get atmosphere model type
    const EphemerisType ephemerisType = getValue< EphemerisType >( settings, keyTree + Keys::type );

    switch ( ephemerisType ) {
    case approximate_planet_positions:
        return boost::make_shared< ApproximatePlanetPositionSettings >(
                    getValue< ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData >(
                        settings, keyTree + Keys::bodyIdentifier ),
                    getValue< bool >( settings, keyTree + Keys::useCircularCoplanarApproximation ) );
    case direct_spice_ephemeris:
    {
        DirectSpiceEphemerisSettings defaults;
        return boost::make_shared< DirectSpiceEphemerisSettings >(
                    getValue( settings, keyTree + Keys::frameOrigin, defaults.getFrameOrigin( ) ),
                    getValue( settings, keyTree + Keys::frameOrientation, defaults.getFrameOrientation( ) ),
                    getValue( settings, keyTree + Keys::correctForStellarAbberation,
                              defaults.getCorrectForStellarAbberation( ) ),
                    getValue( settings, keyTree + Keys::correctForLightTimeAbberation,
                              defaults.getCorrectForLightTimeAbberation( ) ),
                    getValue( settings, keyTree + Keys::convergeLighTimeAbberation,
                              defaults.getConvergeLighTimeAbberation( ) ) );
    }
    case tabulated_ephemeris:
    {
        TabulatedEphemerisSettings defaults( { } );
        TabulatedEphemerisSettings tabulatedEphemerisSettings(
                    getValue< std::map< double, Eigen::Vector6d > >( settings, keyTree + Keys::bodyStateHistory ),
                    getValue( settings, keyTree + Keys::frameOrigin, defaults.getFrameOrigin( ) ),
                    getValue( settings, keyTree + Keys::frameOrientation, defaults.getFrameOrientation( ) ) );
        tabulatedEphemerisSettings.setUseLongDoubleStates(
                    getValue( settings, keyTree + Keys::useLongDoubleStates, defaults.getUseLongDoubleStates( ) ) );
        return boost::make_shared< TabulatedEphemerisSettings >( tabulatedEphemerisSettings );
    }
    case interpolated_spice:
    {
        InterpolatedSpiceEphemerisSettings defaults( TUDAT_NAN, TUDAT_NAN, TUDAT_NAN );
        InterpolatedSpiceEphemerisSettings interpolatedSpiceEphemerisSettings(
                    getEpoch< double >( settings, keyTree + Keys::initialTime ),
                    getEpoch< double >( settings, keyTree + Keys::finalTime ),
                    getNumeric< double >( settings, keyTree + Keys::timeStep ),
                    getValue( settings, keyTree + Keys::frameOrigin, defaults.getFrameOrigin( ) ),
                    getValue( settings, keyTree + Keys::frameOrientation, defaults.getFrameOrientation( ) ),
                    getValue( settings, keyTree + Keys::interpolator, defaults.getInterpolatorSettings( ) ) );
        interpolatedSpiceEphemerisSettings.setUseLongDoubleStates(
                    getValue( settings, keyTree + Keys::useLongDoubleStates, defaults.getUseLongDoubleStates( ) ) );
        return boost::make_shared< InterpolatedSpiceEphemerisSettings >( interpolatedSpiceEphemerisSettings );
    }
    case constant_ephemeris:
    {
        ConstantEphemerisSettings defaults( ( Eigen::Vector6d( ) ) );
        return boost::make_shared< ConstantEphemerisSettings >(
                    getValue< Eigen::Vector6d >( settings, keyTree + Keys::constantState ),
                    getValue( settings, keyTree + Keys::frameOrigin, defaults.getFrameOrigin( ) ),
                    getValue( settings, keyTree + Keys::frameOrientation, defaults.getFrameOrientation( ) ) );
    }
    case kepler_ephemeris:
    {
        KeplerEphemerisSettings defaults( Eigen::Vector6d( ), TUDAT_NAN, TUDAT_NAN );
        return boost::make_shared< KeplerEphemerisSettings >(
                    getValue< Eigen::Vector6d >( settings, keyTree + Keys::initialStateInKeplerianElements ),
                    getEpoch< double >( settings, keyTree + Keys::epochOfInitialState ),
                    getNumeric< double >( settings, keyTree + Keys::centralBodyGravitationalParameter ),
                    getValue( settings, keyTree + Keys::frameOrigin, defaults.getFrameOrigin( ) ),
                    getValue( settings, keyTree + Keys::frameOrientation, defaults.getFrameOrientation( ) ),
                    getValue( settings, keyTree + Keys::rootFinderAbsoluteTolerance,
                              defaults.getRootFinderAbsoluteTolerance( ) ),
                    getValue( settings, keyTree + Keys::rootFinderMaximumNumberOfIterations,
                              defaults.getRootFinderMaximumNumberOfIterations( ) ) );
    }
        /* FIXME
    case custom_ephemeris:
    {

    }
    */
    default:
        throw std::runtime_error( stringFromEnum( ephemerisType, ephemerisTypes )
                                  + " not supported by json_interface." );
    }
}

} // namespace json_interface

} // namespace tudat
