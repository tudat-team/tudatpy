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
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< EphemerisSettings >& ephemerisSettings )
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
        std::shared_ptr< ApproximatePlanetPositionSettings > approximatePlanetPositionSettings =
                std::dynamic_pointer_cast< ApproximatePlanetPositionSettings >( ephemerisSettings );
        assertNonnullptrPointer( approximatePlanetPositionSettings );
        jsonObject[ K::bodyIdentifier ] = approximatePlanetPositionSettings->getBodyIdentifier( );
        jsonObject[ K::useCircularCoplanarApproximation ] =
                approximatePlanetPositionSettings->getUseCircularCoplanarApproximation( );
        return;
    }
    case direct_spice_ephemeris:
    {
        std::shared_ptr< DirectSpiceEphemerisSettings > directSpiceEphemerisSettings =
                std::dynamic_pointer_cast< DirectSpiceEphemerisSettings >( ephemerisSettings );
        assertNonnullptrPointer( directSpiceEphemerisSettings );
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
        std::shared_ptr< InterpolatedSpiceEphemerisSettings > interpolatedSpiceEphemerisSettings =
                std::dynamic_pointer_cast< InterpolatedSpiceEphemerisSettings >( ephemerisSettings );
        assertNonnullptrPointer( interpolatedSpiceEphemerisSettings );
        jsonObject[ K::initialTime ] = interpolatedSpiceEphemerisSettings->getInitialTime( );
        jsonObject[ K::finalTime ] = interpolatedSpiceEphemerisSettings->getFinalTime( );
        jsonObject[ K::timeStep ] = interpolatedSpiceEphemerisSettings->getTimeStep( );
        jsonObject[ K::interpolator ] = interpolatedSpiceEphemerisSettings->getInterpolatorSettings( );
        jsonObject[ K::useLongDoubleStates ] = interpolatedSpiceEphemerisSettings->getUseLongDoubleStates( );
        return;
    }
    case tabulated_ephemeris:
    {
        std::shared_ptr< TabulatedEphemerisSettings > tabulatedEphemerisSettings =
                std::dynamic_pointer_cast< TabulatedEphemerisSettings >( ephemerisSettings );
        assertNonnullptrPointer( tabulatedEphemerisSettings );
        jsonObject[ K::bodyStateHistory ] = tabulatedEphemerisSettings->getBodyStateHistory( );
        jsonObject[ K::useLongDoubleStates ] = tabulatedEphemerisSettings->getUseLongDoubleStates( );
        return;
    }
    case constant_ephemeris:
    {
        std::shared_ptr< ConstantEphemerisSettings > constantEphemerisSettings =
                std::dynamic_pointer_cast< ConstantEphemerisSettings >( ephemerisSettings );
        assertNonnullptrPointer( constantEphemerisSettings );
        jsonObject[ K::constantState ] = constantEphemerisSettings->getConstantState( );
        return;
    }
    case kepler_ephemeris:
    {
        std::shared_ptr< KeplerEphemerisSettings > keplerEphemerisSettings =
                std::dynamic_pointer_cast< KeplerEphemerisSettings >( ephemerisSettings );
        assertNonnullptrPointer( keplerEphemerisSettings );
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
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< EphemerisSettings >& ephemerisSettings )
{
    using namespace ephemerides;
    using namespace json_interface;
    using K = Keys::Body::Ephemeris;

    // Get atmosphere model type
    const EphemerisType ephemerisType = getValue< EphemerisType >( jsonObject, K::type );

    switch ( ephemerisType ) {
    case approximate_planet_positions:
    {
        ephemerisSettings = std::make_shared< ApproximatePlanetPositionSettings >(
                    getValue< ApproximatePlanetPositionsBase::BodiesWithEphemerisData >(
                        jsonObject, K::bodyIdentifier ),
                    getValue< bool >( jsonObject, K::useCircularCoplanarApproximation ) );
        break;
    }
    case direct_spice_ephemeris:
    {
        DirectSpiceEphemerisSettings defaults;
        ephemerisSettings = std::make_shared< DirectSpiceEphemerisSettings >(
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
                    getValue< std::map< double, Eigen::Vector6d > >(
                        jsonObject, K::bodyStateHistory, defaults.getBodyStateHistory( ) ) );
        tabulatedEphemerisSettings.setUseLongDoubleStates(
                    getValue( jsonObject, K::useLongDoubleStates, defaults.getUseLongDoubleStates( ) ) );
        ephemerisSettings = std::make_shared< TabulatedEphemerisSettings >( tabulatedEphemerisSettings );
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
        ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
                    interpolatedSpiceEphemerisSettings );
        break;
    }
    case constant_ephemeris:
    {
        ConstantEphemerisSettings defaults( ( Eigen::Vector6d( ) ) );
        ephemerisSettings = std::make_shared< ConstantEphemerisSettings >(
                    getValue< Eigen::Vector6d >( jsonObject, K::constantState ) );
        break;
    }
    case kepler_ephemeris:
    {
        KeplerEphemerisSettings defaults( Eigen::Vector6d( ), TUDAT_NAN, TUDAT_NAN );
        ephemerisSettings = std::make_shared< KeplerEphemerisSettings >(
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
