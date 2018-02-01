/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "Tudat/JsonInterface/Environment/ephemeris.h"

#include "Tudat/JsonInterface/Mathematics/interpolation.h"

namespace tudat
{

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `EphemerisSettings` object.
void to_json( nlohmann::json& jsonObject, const boost::shared_ptr< EphemerisSettings >& ephemerisSettings )
{
    if ( ! ephemerisSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Body::Ephemeris;

    // Common keys
    const EphemerisType ephemerisType = ephemerisSettings->getEphemerisType( );
    jsonObject[ K::type ] = ephemerisType;
    jsonObject[ K::frameOrigin ] = ephemerisSettings->getFrameOrigin( );
    jsonObject[ K::frameOrientation ] = ephemerisSettings->getFrameOrientation( );
    // jsonObject[ K::makeMultiArc ] = ephemerisSettings->getMakeMultiArcEphemeris( );

    switch ( ephemerisType )
    {
    case approximate_planet_positions:
    {
        boost::shared_ptr< ApproximatePlanetPositionSettings > approximatePlanetPositionSettings =
                boost::dynamic_pointer_cast< ApproximatePlanetPositionSettings >( ephemerisSettings );
        assertNonNullPointer( approximatePlanetPositionSettings );
        jsonObject[ K::bodyIdentifier ] = approximatePlanetPositionSettings->getBodyIdentifier( );
        jsonObject[ K::useCircularCoplanarApproximation ] =
                approximatePlanetPositionSettings->getUseCircularCoplanarApproximation( );
        return;
    }
    case direct_spice_ephemeris:
    {
        boost::shared_ptr< DirectSpiceEphemerisSettings > directSpiceEphemerisSettings =
                boost::dynamic_pointer_cast< DirectSpiceEphemerisSettings >( ephemerisSettings );
        assertNonNullPointer( directSpiceEphemerisSettings );
        jsonObject[ K::correctForStellarAberration ] =
                directSpiceEphemerisSettings->getCorrectForStellarAberration( );
        jsonObject[ K::correctForLightTimeAberration ] =
                directSpiceEphemerisSettings->getCorrectForLightTimeAberration( );
        jsonObject[ K::convergeLighTimeAberration ] =
                directSpiceEphemerisSettings->getConvergeLighTimeAberration( );
        return;
    }
    case interpolated_spice:
    {
        boost::shared_ptr< InterpolatedSpiceEphemerisSettings > interpolatedSpiceEphemerisSettings =
                boost::dynamic_pointer_cast< InterpolatedSpiceEphemerisSettings >( ephemerisSettings );
        assertNonNullPointer( interpolatedSpiceEphemerisSettings );
        jsonObject[ K::initialTime ] = interpolatedSpiceEphemerisSettings->getInitialTime( );
        jsonObject[ K::finalTime ] = interpolatedSpiceEphemerisSettings->getFinalTime( );
        jsonObject[ K::timeStep ] = interpolatedSpiceEphemerisSettings->getTimeStep( );
        jsonObject[ K::interpolator ] = interpolatedSpiceEphemerisSettings->getInterpolatorSettings( );
        jsonObject[ K::useLongDoubleStates ] = interpolatedSpiceEphemerisSettings->getUseLongDoubleStates( );
        return;
    }
    case tabulated_ephemeris:
    {
        boost::shared_ptr< TabulatedEphemerisSettings > tabulatedEphemerisSettings =
                boost::dynamic_pointer_cast< TabulatedEphemerisSettings >( ephemerisSettings );
        assertNonNullPointer( tabulatedEphemerisSettings );
        jsonObject[ K::bodyStateHistory ] = tabulatedEphemerisSettings->getBodyStateHistory( );
        jsonObject[ K::useLongDoubleStates ] = tabulatedEphemerisSettings->getUseLongDoubleStates( );
        return;
    }
    case constant_ephemeris:
    {
        boost::shared_ptr< ConstantEphemerisSettings > constantEphemerisSettings =
                boost::dynamic_pointer_cast< ConstantEphemerisSettings >( ephemerisSettings );
        assertNonNullPointer( constantEphemerisSettings );
        jsonObject[ K::constantState ] = constantEphemerisSettings->getConstantState( );
        return;
    }
    case kepler_ephemeris:
    {
        boost::shared_ptr< KeplerEphemerisSettings > keplerEphemerisSettings =
                boost::dynamic_pointer_cast< KeplerEphemerisSettings >( ephemerisSettings );
        assertNonNullPointer( keplerEphemerisSettings );
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
    default:
        handleUnimplementedEnumValue( ephemerisType, ephemerisTypes, unsupportedEphemerisTypes );
    }
}

//! Create a shared pointer to a `EphemerisSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, boost::shared_ptr< EphemerisSettings >& ephemerisSettings )
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
        break;
    }
    case direct_spice_ephemeris:
    {
        DirectSpiceEphemerisSettings defaults;
        ephemerisSettings = boost::make_shared< DirectSpiceEphemerisSettings >(
                    defaults.getFrameOrigin( ),
                    defaults.getFrameOrientation( ),
                    getValue( jsonObject, K::correctForStellarAberration,
                              defaults.getCorrectForStellarAberration( ) ),
                    getValue( jsonObject, K::correctForLightTimeAberration,
                              defaults.getCorrectForLightTimeAberration( ) ),
                    getValue( jsonObject, K::convergeLighTimeAberration,
                              defaults.getConvergeLighTimeAberration( ) ) );
        break;
    }
    case tabulated_ephemeris:
    {
        TabulatedEphemerisSettings defaults( ( std::map< double, Eigen::Vector6d >( ) ) );
        TabulatedEphemerisSettings tabulatedEphemerisSettings(
                    getValue< std::map< double, Eigen::Vector6d > >( jsonObject, K::bodyStateHistory ) );
        tabulatedEphemerisSettings.setUseLongDoubleStates(
                    getValue( jsonObject, K::useLongDoubleStates, defaults.getUseLongDoubleStates( ) ) );
        ephemerisSettings = boost::make_shared< TabulatedEphemerisSettings >( tabulatedEphemerisSettings );
        break;
    }
    case interpolated_spice:
    {
        InterpolatedSpiceEphemerisSettings defaults( TUDAT_NAN, TUDAT_NAN, TUDAT_NAN );
        InterpolatedSpiceEphemerisSettings interpolatedSpiceEphemerisSettings(
                    getValue< double >( jsonObject, K::initialTime ),
                    getValue< double >( jsonObject, K::finalTime ),
                    getValue< double >( jsonObject, K::timeStep ),
                    defaults.getFrameOrigin( ),
                    defaults.getFrameOrientation( ),
                    getValue( jsonObject, K::interpolator, defaults.getInterpolatorSettings( ) ) );
        interpolatedSpiceEphemerisSettings.setUseLongDoubleStates(
                    getValue( jsonObject, K::useLongDoubleStates, defaults.getUseLongDoubleStates( ) ) );
        ephemerisSettings = boost::make_shared< InterpolatedSpiceEphemerisSettings >(
                    interpolatedSpiceEphemerisSettings );
        break;
    }
    case constant_ephemeris:
    {
        ConstantEphemerisSettings defaults( ( Eigen::Vector6d( ) ) );
        ephemerisSettings = boost::make_shared< ConstantEphemerisSettings >(
                    getValue< Eigen::Vector6d >( jsonObject, K::constantState ) );
        break;
    }
    case kepler_ephemeris:
    {
        KeplerEphemerisSettings defaults( Eigen::Vector6d( ), TUDAT_NAN, TUDAT_NAN );
        ephemerisSettings = boost::make_shared< KeplerEphemerisSettings >(
                    getValue< Eigen::Vector6d >( jsonObject, K::initialStateInKeplerianElements ),
                    getValue< double >( jsonObject, K::epochOfInitialState ),
                    getValue< double >( jsonObject, K::centralBodyGravitationalParameter ),
                    defaults.getFrameOrigin( ),
                    defaults.getFrameOrientation( ),
                    getValue( jsonObject, K::rootFinderAbsoluteTolerance,
                              defaults.getRootFinderAbsoluteTolerance( ) ),
                    getValue( jsonObject, K::rootFinderMaximumNumberOfIterations,
                              defaults.getRootFinderMaximumNumberOfIterations( ) ) );
        break;
    }
    default:
        handleUnimplementedEnumValue( ephemerisType, ephemerisTypes, unsupportedEphemerisTypes );
    }

    if ( isDefined( jsonObject, K::frameOrigin ) )
    {
        ephemerisSettings->resetFrameOrigin( getValue< std::string >( jsonObject, K::frameOrigin ) );
    }

    if ( isDefined( jsonObject, K::frameOrientation ) )
    {
        ephemerisSettings->resetFrameOrientation( getValue< std::string >( jsonObject, K::frameOrientation ) );
    }

    /*
    if ( isDefined( jsonObject, K::makeMultiArc ) )
    {
        ephemerisSettings->resetMakeMultiArcEphemeris( getValue< bool >( jsonObject, K::makeMultiArc ) );
    }
    */
}

} // namespace simulation_setup

} // namespace tudat
