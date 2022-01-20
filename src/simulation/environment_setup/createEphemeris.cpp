/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/lambda/lambda.hpp>

#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/astro/ephemerides/customEphemeris.h"
#include "tudat/astro/ephemerides/keplerEphemeris.h"
#include "tudat/astro/ephemerides/multiArcEphemeris.h"
#include "tudat/astro/ephemerides/tabulatedEphemeris.h"
#include "tudat/astro/ephemerides/approximatePlanetPositions.h"
#include "tudat/astro/ephemerides/approximatePlanetPositionsCircularCoplanar.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/math/interpolators/lagrangeInterpolator.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"

namespace tudat
{

namespace simulation_setup
{

using namespace ephemerides;

//! Function to create an ephemeris model.
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
                std::string inputName = ( directEphemerisSettings->getBodyNameOverride( ) == "" ) ?
                            bodyName : directEphemerisSettings->getBodyNameOverride( );

                // Create corresponding ephemeris object.
                ephemeris = std::make_shared< SpiceEphemeris >(
                            inputName,
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
                std::string inputName = ( interpolatedEphemerisSettings->getBodyNameOverride( ) == "" ) ?
                            bodyName : interpolatedEphemerisSettings->getBodyNameOverride( );
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
#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )

                    ephemeris = createTabulatedEphemerisFromSpice< long double, double >(
                                inputName,
                                static_cast< long double >( interpolatedEphemerisSettings->getInitialTime( ) ),
                                static_cast< long double >( interpolatedEphemerisSettings->getFinalTime( ) ),
                                static_cast< long double >( interpolatedEphemerisSettings->getTimeStep( ) ),
                                interpolatedEphemerisSettings->getFrameOrigin( ),
                                interpolatedEphemerisSettings->getFrameOrientation( ),
                                interpolatedEphemerisSettings->getInterpolatorSettings( ) );
#else
                    throw std::runtime_error( "Error, long double compilation is turned off; requested long double tabulated ephemeris" );
#endif
                }
            }
            break;
        }
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
#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )

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
                    throw std::runtime_error( "Error, long double compilation is turned off; requested long double tabulated ephemeris" );
#endif
                }
            }
            break;
        }
        case auto_generated_tabulated_ephemeris:
        {
            // Check consistency of type and class.
            std::shared_ptr< AutoGeneratedTabulatedEphemerisSettings > tabulatedEphemerisSettings =
                    std::dynamic_pointer_cast< AutoGeneratedTabulatedEphemerisSettings >( ephemerisSettings );
            if( tabulatedEphemerisSettings == nullptr )
            {
                throw std::runtime_error( "Error, expected auto-generate tabulated ephemeris settings for " + bodyName );
            }
            else
            {
                // Create ephemeris
                ephemeris = getTabulatedEphemeris(
                            createBodyEphemeris( tabulatedEphemerisSettings->getEphemerisSettings( ), bodyName ),
                            tabulatedEphemerisSettings->getStartTime( ), tabulatedEphemerisSettings->getEndTime( ),
                            tabulatedEphemerisSettings->getTimeStep( ), tabulatedEphemerisSettings->getInterpolatorSettings( ) );
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
                throw std::runtime_error( "Error, expected custom ephemeris settings for " + bodyName );
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
            std::shared_ptr< ApproximateJplEphemerisSettings > approximateEphemerisSettings =
                    std::dynamic_pointer_cast< ApproximateJplEphemerisSettings >(
                        ephemerisSettings );
            if( approximateEphemerisSettings == nullptr )
            {
                throw std::runtime_error(
                            "Error, expected approximate ephemeris settings for body " + bodyName );
            }
            else
            {

//                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData bodyToUse;
//                if( approximateEphemerisSettings->getBodyIdentifier( ) ==
//                        ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::undefined )
//                {
//                    try
//                    {
//                        bodyToUse = ephemerides::ApproximatePlanetPositionsBase::getBodiesWithEphemerisDataId( bodyName );
//                    }
//                    catch( std::runtime_error const& )

//                    {
//                        throw std::runtime_error( "Error, approximate ephemeris not available for body: " + bodyName + " when creating ephemeris." );
//                    }
//                }
//                else
//                {
//                    bodyToUse = approximateEphemerisSettings->getBodyIdentifier( );
//                }
//=======
//>>>>>>> origin/feature/mga_estimation_refactor_merge

                // Create corresponding ephemeris object.
                if( approximateEphemerisSettings->getUseCircularCoplanarApproximation( ) )
                {
                    ephemeris = std::make_shared< ApproximateJplCircularCoplanarEphemeris >(
                                approximateEphemerisSettings->getBodyName( ) );
                }
                else
                {
                    ephemeris = std::make_shared< ApproximateJplEphemeris >(
                                approximateEphemerisSettings->getBodyName( ) );
                }
            }
            break;
        }
		case direct_tle_ephemeris:
		{
			// Check consistency of type and class.
			std::shared_ptr< DirectTleEphemerisSettings > directTleEphemerisSettings =
					std::dynamic_pointer_cast< DirectTleEphemerisSettings >( ephemerisSettings );
			if( directTleEphemerisSettings == nullptr )
			{
				throw std::runtime_error(
						"Error, expected direct TLE ephemeris settings for body " + bodyName );
			}
			// Check if the Earth is present in the simulation
			else
			{
				//std::string inputName = bodyName;

				// Check period of the satellite for correct SDP setting


				// Create corresponding ephemeris object.
				ephemeris = std::make_shared< TleEphemeris >(
						directTleEphemerisSettings->getFrameOrigin(),
						directTleEphemerisSettings->getFrameOrientation(),
						directTleEphemerisSettings->getTle()
						);
			}
			break;
		}
		case interpolated_tle_ephemeris:
		{
			// Check consistency of type and class.
			std::shared_ptr< InterpolatedTleEphemerisSettings > interpolatedTleEphemerisSettings =
					std::dynamic_pointer_cast< InterpolatedTleEphemerisSettings >(
							ephemerisSettings );
			if( interpolatedTleEphemerisSettings == nullptr )
			{
				throw std::runtime_error(
						"Error, expected interpolated TLE ephemeris settings for body " + bodyName );
			}
			else
			{

				ephemeris = createTabulatedEphemerisFromTLE< double, double >(
						bodyName,
						interpolatedTleEphemerisSettings->getInitialTime( ),
						interpolatedTleEphemerisSettings->getFinalTime( ),
						interpolatedTleEphemerisSettings->getTimeStep( ),
						interpolatedTleEphemerisSettings->getFrameOrigin( ),
						interpolatedTleEphemerisSettings->getFrameOrientation( ),
						interpolatedTleEphemerisSettings->getTle( ),
						interpolatedTleEphemerisSettings->getInterpolatorSettings( )
						);
			}
			break;
		}
        case scaled_ephemeris:
        {
            // Check consistency of type and class.
            std::shared_ptr< ScaledEphemerisSettings > scaledEphemeriSettings =
                    std::dynamic_pointer_cast< ScaledEphemerisSettings >(
                        ephemerisSettings );
            if( scaledEphemeriSettings == nullptr )
            {
                throw std::runtime_error(
                            "Error, expected scaled ephemeris settings for body " + bodyName );
            }
            else
            {
                std::shared_ptr< Ephemeris > baseEphemeris = createBodyEphemeris(
                            scaledEphemeriSettings->getBaseSettings( ), bodyName );
                ephemeris = std::make_shared< ScaledEphemeris >(
                            baseEphemeris, scaledEphemeriSettings->getScaling( ), scaledEphemeriSettings->getIsScalingAbsolute( ) );
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
