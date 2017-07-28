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
    jsonObject = json_interface::stringFromEnum( ephemerisType, ephemerisTypes );
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
        using K = Keys::Body::Ephemeris;

        // Common keys
        jsonObject[ K::type ] = ephemerisSettings->getEphemerisType( );
        jsonObject[ K::frameOrigin ] = ephemerisSettings->getFrameOrigin( );
        jsonObject[ K::frameOrientation ] = ephemerisSettings->getFrameOrientation( );
        jsonObject[ K::makeMultiArc ] = ephemerisSettings->getMakeMultiArcEphemeris( );

        /// ApproximatePlanetPositionSettings
        boost::shared_ptr< ApproximatePlanetPositionSettings > approximatePlanetPositionSettings =
                boost::dynamic_pointer_cast< ApproximatePlanetPositionSettings >( ephemerisSettings );
        if ( approximatePlanetPositionSettings )
        {
            jsonObject[ K::bodyIdentifier ] = approximatePlanetPositionSettings->getBodyIdentifier( );
            jsonObject[ K::useCircularCoplanarApproximation ] =
                    approximatePlanetPositionSettings->getUseCircularCoplanarApproximation( );
            return;
        }

        /// DirectSpiceEphemerisSettings
        boost::shared_ptr< DirectSpiceEphemerisSettings > directSpiceEphemerisSettings =
                boost::dynamic_pointer_cast< DirectSpiceEphemerisSettings >( ephemerisSettings );
        if ( directSpiceEphemerisSettings )
        {
            jsonObject[ K::correctForStellarAbberation ] =
                    directSpiceEphemerisSettings->getCorrectForStellarAbberation( );
            jsonObject[ K::correctForLightTimeAbberation ] =
                    directSpiceEphemerisSettings->getCorrectForLightTimeAbberation( );
            jsonObject[ K::convergeLighTimeAbberation ] =
                    directSpiceEphemerisSettings->getConvergeLighTimeAbberation( );

            /// InterpolatedSpiceEphemerisSettings
            boost::shared_ptr< InterpolatedSpiceEphemerisSettings > interpolatedSpiceEphemerisSettings =
                    boost::dynamic_pointer_cast< InterpolatedSpiceEphemerisSettings >( directSpiceEphemerisSettings );
            if ( interpolatedSpiceEphemerisSettings )
            {
                jsonObject[ K::initialTime ] = interpolatedSpiceEphemerisSettings->getInitialTime( );
                jsonObject[ K::finalTime ] = interpolatedSpiceEphemerisSettings->getFinalTime( );
                jsonObject[ K::timeStep ] = interpolatedSpiceEphemerisSettings->getTimeStep( );
                jsonObject[ K::interpolator ] = interpolatedSpiceEphemerisSettings->getInterpolatorSettings( );
                jsonObject[ K::useLongDoubleStates ] = interpolatedSpiceEphemerisSettings->getUseLongDoubleStates( );
                return;
            }

            return;
        }

        /// TabulatedEphemerisSettings
        boost::shared_ptr< TabulatedEphemerisSettings > tabulatedEphemerisSettings =
                boost::dynamic_pointer_cast< TabulatedEphemerisSettings >( ephemerisSettings );
        if ( tabulatedEphemerisSettings )
        {
            jsonObject[ K::bodyStateHistory ] = tabulatedEphemerisSettings->getBodyStateHistory( );
            jsonObject[ K::useLongDoubleStates ] = tabulatedEphemerisSettings->getUseLongDoubleStates( );
            return;
        }

        /// ConstantEphemerisSettings
        boost::shared_ptr< ConstantEphemerisSettings > constantEphemerisSettings =
                boost::dynamic_pointer_cast< ConstantEphemerisSettings >( ephemerisSettings );
        if ( constantEphemerisSettings )
        {
            jsonObject[ K::constantState ] = constantEphemerisSettings->getConstantState( );
            return;
        }

        /// KeplerEphemerisSettings
        boost::shared_ptr< KeplerEphemerisSettings > keplerEphemerisSettings =
                boost::dynamic_pointer_cast< KeplerEphemerisSettings >( ephemerisSettings );
        if ( keplerEphemerisSettings )
        {
            jsonObject[ K::initialStateInKeplerianElements ] =
                    keplerEphemerisSettings->getInitialStateInKeplerianElements( );
            jsonObject[ K::epochOfInitialState ] =
                    keplerEphemerisSettings->getEpochOfInitialState( );
            jsonObject[ K::centralBodyGravitationalParameter ] =
                    keplerEphemerisSettings->getCentralBodyGravitationalParameter( );
            jsonObject[ K::rootFinderAbsoluteTolerance ] =
                    keplerEphemerisSettings->getRootFinderAbsoluteTolerance( );
            jsonObject[ K::rootFinderMaximumNumberOfIterations ] =
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

//! Create a shared pointer to a `EphemerisSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< EphemerisSettings >& ephemerisSettings )
{
    using namespace ephemerides;
    using namespace json_interface;
    using K = Keys::Body::Ephemeris;

    // Get atmosphere model type
    const EphemerisType ephemerisType = getValue< EphemerisType >( jsonObject, K::type );

    switch ( ephemerisType ) {
    case approximate_planet_positions:
    {
        ephemerisSettings = boost::make_shared< ApproximatePlanetPositionSettings >(
                    getValue< ApproximatePlanetPositionsBase::BodiesWithEphemerisData >(
                        jsonObject, K::bodyIdentifier ),
                    getValue< bool >( jsonObject, K::useCircularCoplanarApproximation ) );
        return;
    }
    case direct_spice_ephemeris:
    {
        DirectSpiceEphemerisSettings defaults;
        ephemerisSettings = boost::make_shared< DirectSpiceEphemerisSettings >(
                    getValue( jsonObject, K::frameOrigin, defaults.getFrameOrigin( ) ),
                    getValue( jsonObject, K::frameOrientation, defaults.getFrameOrientation( ) ),
                    getValue( jsonObject, K::correctForStellarAbberation,
                              defaults.getCorrectForStellarAbberation( ) ),
                    getValue( jsonObject, K::correctForLightTimeAbberation,
                              defaults.getCorrectForLightTimeAbberation( ) ),
                    getValue( jsonObject, K::convergeLighTimeAbberation,
                              defaults.getConvergeLighTimeAbberation( ) ) );
        return;
    }
    case tabulated_ephemeris:
    {
        /* FIXME::MAP
        TabulatedEphemerisSettings defaults( { } );
        TabulatedEphemerisSettings tabulatedEphemerisSettings(
                    getValue< std::map< double, Eigen::Vector6d > >( jsonObject, K::bodyStateHistory ),
                    getValue( jsonObject, K::frameOrigin, defaults.getFrameOrigin( ) ),
                    getValue( jsonObject, K::frameOrientation, defaults.getFrameOrientation( ) ) );
        tabulatedEphemerisSettings.setUseLongDoubleStates(
                    getValue( jsonObject, K::useLongDoubleStates, defaults.getUseLongDoubleStates( ) ) );
        ephemerisSettings = boost::make_shared< TabulatedEphemerisSettings >( tabulatedEphemerisSettings );
        */
        return;
    }
    case interpolated_spice:
    {
        InterpolatedSpiceEphemerisSettings defaults( TUDAT_NAN, TUDAT_NAN, TUDAT_NAN );
        InterpolatedSpiceEphemerisSettings interpolatedSpiceEphemerisSettings(
                    getEpoch< double >( jsonObject, K::initialTime ),
                    getEpoch< double >( jsonObject, K::finalTime ),
                    getNumeric< double >( jsonObject, K::timeStep ),
                    getValue( jsonObject, K::frameOrigin, defaults.getFrameOrigin( ) ),
                    getValue( jsonObject, K::frameOrientation, defaults.getFrameOrientation( ) ),
                    getValue( jsonObject, K::interpolator, defaults.getInterpolatorSettings( ) ) );
        interpolatedSpiceEphemerisSettings.setUseLongDoubleStates(
                    getValue( jsonObject, K::useLongDoubleStates, defaults.getUseLongDoubleStates( ) ) );
        ephemerisSettings = boost::make_shared< InterpolatedSpiceEphemerisSettings >(
                    interpolatedSpiceEphemerisSettings );
        return;
    }
    case constant_ephemeris:
    {
        ConstantEphemerisSettings defaults( ( Eigen::Vector6d( ) ) );
        ephemerisSettings = boost::make_shared< ConstantEphemerisSettings >(
                    getValue< Eigen::Vector6d >( jsonObject, K::constantState ),
                    getValue( jsonObject, K::frameOrigin, defaults.getFrameOrigin( ) ),
                    getValue( jsonObject, K::frameOrientation, defaults.getFrameOrientation( ) ) );
        return;
    }
    case kepler_ephemeris:
    {
        KeplerEphemerisSettings defaults( Eigen::Vector6d( ), TUDAT_NAN, TUDAT_NAN );
        ephemerisSettings = boost::make_shared< KeplerEphemerisSettings >(
                    getValue< Eigen::Vector6d >( jsonObject, K::initialStateInKeplerianElements ),
                    getEpoch< double >( jsonObject, K::epochOfInitialState ),
                    getNumeric< double >( jsonObject, K::centralBodyGravitationalParameter ),
                    getValue( jsonObject, K::frameOrigin, defaults.getFrameOrigin( ) ),
                    getValue( jsonObject, K::frameOrientation, defaults.getFrameOrientation( ) ),
                    getValue( jsonObject, K::rootFinderAbsoluteTolerance,
                              defaults.getRootFinderAbsoluteTolerance( ) ),
                    getValue( jsonObject, K::rootFinderMaximumNumberOfIterations,
                              defaults.getRootFinderMaximumNumberOfIterations( ) ) );
        return;
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

} // namespace simulation_setup

} // namespace tudat
