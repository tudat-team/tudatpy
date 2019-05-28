/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/lambda/lambda.hpp>
#if USE_CSPICE
#include "Tudat/External/SpiceInterface/spiceEphemeris.h"
#endif

#include "Tudat/Astrodynamics/Ephemerides/customEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/keplerEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/multiArcEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/tabulatedEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositionsCircularCoplanar.h"
#include "Tudat/Astrodynamics/Ephemerides/constantEphemeris.h"
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.h"

namespace tudat
{

namespace simulation_setup
{

using namespace ephemerides;

//! Function to create a ephemeris model.
std::shared_ptr< ephemerides::Ephemeris > createBodyEphemeris(
        const std::shared_ptr< EphemerisSettings > ephemerisSettings,
        const std::string& bodyName )
{
    // Declare return object.
    std::shared_ptr< ephemerides::Ephemeris > ephemeris;

    if( ephemerisSettings->getMakeMultiArcEphemeris( ) )
    {
        std::map< double, std::shared_ptr< Ephemeris > > singleArcEphemerides;
        ephemerisSettings->resetMakeMultiArcEphemeris( false );

        singleArcEphemerides[ -std::numeric_limits< double >::lowest( ) ] = createBodyEphemeris(
                    ephemerisSettings, bodyName );

        ephemeris = std::make_shared< MultiArcEphemeris >(
                    singleArcEphemerides, ephemerisSettings->getFrameOrigin( ), ephemerisSettings->getFrameOrientation( ) );
    }
    else
    {

        // Check which type of ephemeris model is to be created.
        switch( ephemerisSettings->getEphemerisType( ) )
        {
#if USE_CSPICE
        case direct_spice_ephemeris:
        {
            // Check consistency of type and class.
            std::shared_ptr< DirectSpiceEphemerisSettings > directEphemerisSettings =
                    std::dynamic_pointer_cast< DirectSpiceEphemerisSettings >( ephemerisSettings );
            if( directEphemerisSettings == nullptr )
            {
                throw std::runtime_error(
                            "Error, expected direct spice ephemeris settings for body " + bodyName );
            }
            else
            {
                // Create corresponding ephemeris object.
                ephemeris = std::make_shared< SpiceEphemeris >(
                            bodyName,
                            directEphemerisSettings->getFrameOrigin( ),
                            directEphemerisSettings->getCorrectForStellarAberration( ),
                            directEphemerisSettings->getCorrectForLightTimeAberration( ),
                            directEphemerisSettings->getConvergeLighTimeAberration( ),
                            directEphemerisSettings->getFrameOrientation( ) );
            }
            break;
        }
        case interpolated_spice:
        {
            // Check consistency of type and class.
            std::shared_ptr< InterpolatedSpiceEphemerisSettings > interpolatedEphemerisSettings =
                    std::dynamic_pointer_cast< InterpolatedSpiceEphemerisSettings >(
                        ephemerisSettings );
            if( interpolatedEphemerisSettings == nullptr )
            {
                throw std::runtime_error(
                            "Error, expected tabulated spice ephemeris settings for body " + bodyName );
            }
            else
            {
                // Since only the barycenters of planetary systems are included in the standard DE
                // ephemerides, append 'Barycenter' to body name.
                std::string inputName;
                inputName = bodyName;
                if( bodyName == "Mars" ||
                        bodyName == "Jupiter"  || bodyName == "Saturn" ||
                        bodyName == "Uranus" || bodyName == "Neptune" )
                {
                    inputName += " Barycenter";
                    std::cerr << "Warning, position of " << bodyName << " taken as barycenter of that body's "
                              << "planetary system." << std::endl;
                }

                // Create corresponding ephemeris object.
                if( !interpolatedEphemerisSettings->getUseLongDoubleStates( ) )
                {
                    ephemeris = createTabulatedEphemerisFromSpice< double, double >(
                                inputName,
                                interpolatedEphemerisSettings->getInitialTime( ),
                                interpolatedEphemerisSettings->getFinalTime( ),
                                interpolatedEphemerisSettings->getTimeStep( ),
                                interpolatedEphemerisSettings->getFrameOrigin( ),
                                interpolatedEphemerisSettings->getFrameOrientation( ),
                                interpolatedEphemerisSettings->getInterpolatorSettings( ) );
                }
                else
                {
#if( BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS )

                    ephemeris = createTabulatedEphemerisFromSpice< long double, double >(
                                inputName,
                                static_cast< long double >( interpolatedEphemerisSettings->getInitialTime( ) ),
                                static_cast< long double >( interpolatedEphemerisSettings->getFinalTime( ) ),
                                static_cast< long double >( interpolatedEphemerisSettings->getTimeStep( ) ),
                                interpolatedEphemerisSettings->getFrameOrigin( ),
                                interpolatedEphemerisSettings->getFrameOrientation( ),
                                interpolatedEphemerisSettings->getInterpolatorSettings( ) );
#else
                    throw std::runtime_error( "Error, long double compilation is turned off; requested long doubel tabulated ephemeris" );
#endif
                }
            }
            break;
        }
#endif
        case tabulated_ephemeris:
        {
            // Check consistency of type and class.
            std::shared_ptr< TabulatedEphemerisSettings > tabulatedEphemerisSettings =
                    std::dynamic_pointer_cast< TabulatedEphemerisSettings >( ephemerisSettings );
            if( tabulatedEphemerisSettings == nullptr )
            {
                throw std::runtime_error(
                            "Error, expected tabulated ephemeris settings for body " + bodyName );
            }
            else
            {
                // Create corresponding ephemeris object.
                if( !tabulatedEphemerisSettings->getUseLongDoubleStates( ) )
                {
                    if( tabulatedEphemerisSettings->getBodyStateHistory( ).size( ) != 0 )
                    {
                        ephemeris = std::make_shared< TabulatedCartesianEphemeris< > >(
                                    std::make_shared<
                                    interpolators::LagrangeInterpolator< double, Eigen::Vector6d > >
                                    ( tabulatedEphemerisSettings->getBodyStateHistory( ), 6,
                                      interpolators::huntingAlgorithm,
                                      interpolators::lagrange_cubic_spline_boundary_interpolation ),
                                    tabulatedEphemerisSettings->getFrameOrigin( ),
                                    tabulatedEphemerisSettings->getFrameOrientation( ) );
                    }
                    else
                    {
                        ephemeris = std::make_shared< TabulatedCartesianEphemeris< > >(
                                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d > >( ),
                                      tabulatedEphemerisSettings->getFrameOrigin( ),
                                      tabulatedEphemerisSettings->getFrameOrientation( ) );
                    }
                }
                else
                {
#if( BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS )

                    // Cast input history to required type.
                    if( tabulatedEphemerisSettings->getBodyStateHistory( ).size( ) != 0 )
                    {
                        std::map< double, Eigen::Vector6d > originalStateHistory =
                                tabulatedEphemerisSettings->getBodyStateHistory( );
                        std::map< double, Eigen::Matrix< long double, 6, 1 > > longStateHistory;

                        for( std::map< double, Eigen::Vector6d >::const_iterator stateIterator =
                             originalStateHistory.begin( ); stateIterator != originalStateHistory.end( ); stateIterator++ )
                        {
                            longStateHistory[ stateIterator->first ] = stateIterator->second.cast< long double >( );
                            ephemeris =
                                    std::make_shared< TabulatedCartesianEphemeris< long double, double > >(
                                        std::make_shared< interpolators::LagrangeInterpolator<
                                        double, Eigen::Matrix< long double, 6, 1 > > >
                                        ( longStateHistory, 6,
                                          interpolators::huntingAlgorithm,
                                          interpolators::lagrange_cubic_spline_boundary_interpolation ),
                                        tabulatedEphemerisSettings->getFrameOrigin( ),
                                        tabulatedEphemerisSettings->getFrameOrientation( ) );
                        }
                    }
                    else
                    {
                        ephemeris = std::make_shared< TabulatedCartesianEphemeris< long double, double > >(
                                    std::shared_ptr< interpolators::OneDimensionalInterpolator<
                                    double, Eigen::Matrix< long double, 6, 1 > > >( ),
                                      tabulatedEphemerisSettings->getFrameOrigin( ),
                                      tabulatedEphemerisSettings->getFrameOrientation( ) );
                    }
#else
                    throw std::runtime_error( "Error, long double compilation is turned off; requested long doubel tabulated ephemeris" );
#endif
                }
            }
            break;
        }
        case constant_ephemeris:
        {
            // Check consistency of type and class.
            std::shared_ptr< ConstantEphemerisSettings > constantEphemerisSettings =
                    std::dynamic_pointer_cast< ConstantEphemerisSettings >( ephemerisSettings );
            if( constantEphemerisSettings == nullptr )
            {
                throw std::runtime_error( "Error, expected constant ephemeris settings for " + bodyName );
            }
            else
            {
                // Create ephemeris
                ephemeris = std::make_shared< ConstantEphemeris >(
                            [ = ]( ){ return constantEphemerisSettings->getConstantState( ); },
                            constantEphemerisSettings->getFrameOrigin( ),
                            constantEphemerisSettings->getFrameOrientation( ) );
            }
            break;
        }
        case custom_ephemeris:
        {
            // Check consistency of type and class.
            std::shared_ptr< CustomEphemerisSettings > customEphemerisSettings =
                    std::dynamic_pointer_cast< CustomEphemerisSettings >( ephemerisSettings );
            if( customEphemerisSettings == nullptr )
            {
                throw std::runtime_error( "Error, expected constant ephemeris settings for " + bodyName );
            }
            else
            {
                // Create ephemeris
                ephemeris = std::make_shared< CustomEphemeris >(
                            customEphemerisSettings->getCustomStateFunction( ),
                            customEphemerisSettings->getFrameOrigin( ),
                            customEphemerisSettings->getFrameOrientation( ) );
            }
            break;
        }

        case kepler_ephemeris:
        {
            // Check consistency of type and class.
            std::shared_ptr< KeplerEphemerisSettings > keplerEphemerisSettings =
                    std::dynamic_pointer_cast< KeplerEphemerisSettings >( ephemerisSettings );
            if( keplerEphemerisSettings == nullptr )
            {
                throw std::runtime_error( "Error, expected Kepler ephemeris settings for " + bodyName );
            }
            else
            {
                // Create ephemeris
                ephemeris = std::make_shared< KeplerEphemeris >(
                            keplerEphemerisSettings->getInitialStateInKeplerianElements( ),
                            keplerEphemerisSettings->getEpochOfInitialState( ),
                            keplerEphemerisSettings->getCentralBodyGravitationalParameter( ),
                            keplerEphemerisSettings->getFrameOrigin( ),
                            keplerEphemerisSettings->getFrameOrientation( ),
                            keplerEphemerisSettings->getRootFinderAbsoluteTolerance( ),
                            keplerEphemerisSettings->getRootFinderMaximumNumberOfIterations( ) );
            }
            break;
        }
        case approximate_planet_positions:
        {
            // Check consistency of type and class.
            std::shared_ptr< ApproximatePlanetPositionSettings > approximateEphemerisSettings =
                    std::dynamic_pointer_cast< ApproximatePlanetPositionSettings >(
                        ephemerisSettings );
            if( approximateEphemerisSettings == nullptr )
            {
                throw std::runtime_error(
                            "Error, expected approximate ephemeris settings for body " + bodyName );
            }
            else
            {
                // Create corresponding ephemeris object.
                if( approximateEphemerisSettings->getUseCircularCoplanarApproximation( ) )
                {
                    ephemeris = std::make_shared< ApproximatePlanetPositionsCircularCoplanar >(
                                approximateEphemerisSettings->getBodyIdentifier( ) );
                }
                else
                {
                    ephemeris = std::make_shared< ApproximatePlanetPositions >(
                                approximateEphemerisSettings->getBodyIdentifier( ) );
                }
            }
            break;
        }
        default:
        {
            throw std::runtime_error(
                        "Error, did not recognize ephemeris model settings type " +
                        std::to_string( ephemerisSettings->getEphemerisType( ) ) );
        }
        }
    }
    return ephemeris;

}

//! Function that retrieves the time interval at which an ephemeris can be safely interrogated
std::pair< double, double > getSafeInterpolationInterval( const std::shared_ptr< ephemerides::Ephemeris > ephemerisModel )
{
    // Make default output pair
    std::pair< double, double > safeInterval = std::make_pair(
                std::numeric_limits< double >::lowest( ),  std::numeric_limits< double >::max( ) );

    // Check if model is tabulated, and retrieve safe interval from model
    if( isTabulatedEphemeris( ephemerisModel ) )
    {
        safeInterval = getTabulatedEphemerisSafeInterval( ephemerisModel );
    }
    // Check if model is multi-arc, and retrieve safe intervals from first and last arc.
    else if( std::dynamic_pointer_cast< ephemerides::MultiArcEphemeris >( ephemerisModel ) != nullptr )
    {
        std::shared_ptr< ephemerides::MultiArcEphemeris > multiArcEphemerisModel  =
                std::dynamic_pointer_cast< ephemerides::MultiArcEphemeris >( ephemerisModel );
        safeInterval.first = getSafeInterpolationInterval( multiArcEphemerisModel->getSingleArcEphemerides( ).at( 0 ) ).first;
        safeInterval.second = getSafeInterpolationInterval(
                    multiArcEphemerisModel->getSingleArcEphemerides( ).at(
                        multiArcEphemerisModel->getSingleArcEphemerides( ).size( ) - 1 ) ).second;
    }
    return safeInterval;
}


} // namespace simulation_setup

} // namespace tudat
