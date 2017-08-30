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

#include "Tudat/External/JsonInterface/Mathematics/interpolation.h"

namespace tudat
{

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `EphemerisSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< EphemerisSettings >& ephemerisSettings )
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
    jsonObject[ K::makeMultiArc ] = ephemerisSettings->getMakeMultiArcEphemeris( );

    switch ( ephemerisType )
    {
    case approximate_planet_positions:
    {
        boost::shared_ptr< ApproximatePlanetPositionSettings > approximatePlanetPositionSettings =
                boost::dynamic_pointer_cast< ApproximatePlanetPositionSettings >( ephemerisSettings );
        enforceNonNullPointer( approximatePlanetPositionSettings );
        jsonObject[ K::bodyIdentifier ] = approximatePlanetPositionSettings->getBodyIdentifier( );
        jsonObject[ K::useCircularCoplanarApproximation ] =
                approximatePlanetPositionSettings->getUseCircularCoplanarApproximation( );
        return;
    }
    case direct_spice_ephemeris:
    case interpolated_spice:
    {
        boost::shared_ptr< DirectSpiceEphemerisSettings > directSpiceEphemerisSettings =
                boost::dynamic_pointer_cast< DirectSpiceEphemerisSettings >( ephemerisSettings );
        enforceNonNullPointer( directSpiceEphemerisSettings );
        jsonObject[ K::correctForStellarAbberation ] =
                directSpiceEphemerisSettings->getCorrectForStellarAbberation( );
        jsonObject[ K::correctForLightTimeAbberation ] =
                directSpiceEphemerisSettings->getCorrectForLightTimeAbberation( );
        jsonObject[ K::convergeLighTimeAbberation ] =
                directSpiceEphemerisSettings->getConvergeLighTimeAbberation( );

        if ( ephemerisType == interpolated_spice )
        {
            boost::shared_ptr< InterpolatedSpiceEphemerisSettings > interpolatedSpiceEphemerisSettings =
                    boost::dynamic_pointer_cast< InterpolatedSpiceEphemerisSettings >( ephemerisSettings );
            enforceNonNullPointer( interpolatedSpiceEphemerisSettings );
            jsonObject[ K::initialTime ] = interpolatedSpiceEphemerisSettings->getInitialTime( );
            jsonObject[ K::finalTime ] = interpolatedSpiceEphemerisSettings->getFinalTime( );
            jsonObject[ K::timeStep ] = interpolatedSpiceEphemerisSettings->getTimeStep( );
            jsonObject[ K::interpolator ] = interpolatedSpiceEphemerisSettings->getInterpolatorSettings( );
            jsonObject[ K::useLongDoubleStates ] = interpolatedSpiceEphemerisSettings->getUseLongDoubleStates( );
            return;
        }

        return;
    }
    case tabulated_ephemeris:
    {
        boost::shared_ptr< TabulatedEphemerisSettings > tabulatedEphemerisSettings =
                boost::dynamic_pointer_cast< TabulatedEphemerisSettings >( ephemerisSettings );
        enforceNonNullPointer( tabulatedEphemerisSettings );
        jsonObject[ K::bodyStateHistory ] = tabulatedEphemerisSettings->getBodyStateHistory( );
        jsonObject[ K::useLongDoubleStates ] = tabulatedEphemerisSettings->getUseLongDoubleStates( );
        return;
    }
    case constant_ephemeris:
    {
        boost::shared_ptr< ConstantEphemerisSettings > constantEphemerisSettings =
                boost::dynamic_pointer_cast< ConstantEphemerisSettings >( ephemerisSettings );
        enforceNonNullPointer( constantEphemerisSettings );
        jsonObject[ K::constantState ] = constantEphemerisSettings->getConstantState( );
        return;
    }
    case kepler_ephemeris:
    {
        boost::shared_ptr< KeplerEphemerisSettings > keplerEphemerisSettings =
                boost::dynamic_pointer_cast< KeplerEphemerisSettings >( ephemerisSettings );
        enforceNonNullPointer( keplerEphemerisSettings );
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
        break;
    }
    case direct_spice_ephemeris:
    {
        DirectSpiceEphemerisSettings defaults;
        ephemerisSettings = boost::make_shared< DirectSpiceEphemerisSettings >(
                    defaults.getFrameOrigin( ),
                    defaults.getFrameOrientation( ),
                    getValue( jsonObject, K::correctForStellarAbberation,
                              defaults.getCorrectForStellarAbberation( ) ),
                    getValue( jsonObject, K::correctForLightTimeAbberation,
                              defaults.getCorrectForLightTimeAbberation( ) ),
                    getValue( jsonObject, K::convergeLighTimeAbberation,
                              defaults.getConvergeLighTimeAbberation( ) ) );
        break;
    }
    case tabulated_ephemeris:
    {
        TabulatedEphemerisSettings defaults( { } );
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
                    getEpoch< double >( jsonObject, K::initialTime ),
                    getEpoch< double >( jsonObject, K::finalTime ),
                    getNumeric< double >( jsonObject, K::timeStep ),
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
                    getEpoch< double >( jsonObject, K::epochOfInitialState ),
                    getNumeric< double >( jsonObject, K::centralBodyGravitationalParameter ),
                    defaults.getFrameOrigin( ),
                    defaults.getFrameOrientation( ),
                    getValue( jsonObject, K::rootFinderAbsoluteTolerance,
                              defaults.getRootFinderAbsoluteTolerance( ) ),
                    getValue( jsonObject, K::rootFinderMaximumNumberOfIterations,
                              defaults.getRootFinderMaximumNumberOfIterations( ) ) );
        break;
    }
        /* FIXME
    case custom_ephemeris:
    {

    }
    */
    default:
        handleUnimplementedEnumValue( ephemerisType, ephemerisTypes, unsupportedEphemerisTypes );
    }

    if ( defined( jsonObject, K::frameOrigin ) )
    {
        ephemerisSettings->resetFrameOrigin( getValue< std::string >( jsonObject, K::frameOrigin ) );
    }

    if ( defined( jsonObject, K::frameOrientation ) )
    {
        ephemerisSettings->resetFrameOrientation( getValue< std::string >( jsonObject, K::frameOrientation ) );
    }

    if ( defined( jsonObject, K::makeMultiArc ) )
    {
        ephemerisSettings->resetMakeMultiArcEphemeris( getValue< bool >( jsonObject, K::makeMultiArc ) );
    }
}

} // namespace simulation_setup

} // namespace tudat
