/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/lambda/lambda.hpp>
#include <boost/lexical_cast.hpp>

#if USE_CSPICE
#include "Tudat/External/SpiceInterface/spiceEphemeris.h"
#endif

#include "Tudat/Astrodynamics/Ephemerides/keplerEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/tabulatedEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositionsCircularCoplanar.h"
#include "Tudat/Astrodynamics/Ephemerides/constantEphemeris.h"
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"
#include "Tudat/SimulationSetup/createEphemeris.h"

namespace tudat
{

namespace simulation_setup
{

using namespace ephemerides;

#if USE_CSPICE
//! Function to create a tabulated ephemeris using data from Spice.
boost::shared_ptr< Ephemeris > createTabulatedEphemerisFromSpice(
        const std::string& body,
        const double initialTime,
        const double endTime,
        const double timeStep,
        const std::string& observerName,
        const std::string& referenceFrameName,
        boost::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
        boost::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ) )
{
    using namespace interpolators;

    std::map< double, basic_mathematics::Vector6d > timeHistoryOfState;

    // Calculate state from spice at given time intervals and store in timeHistoryOfState.
    double currentTime = initialTime;
    while( currentTime < endTime )
    {
        timeHistoryOfState[ currentTime ] = spice_interface::getBodyCartesianStateAtEpoch(
                    body, observerName, referenceFrameName, "none", currentTime );
        currentTime += timeStep;
    }

    // Create interpolator.
    boost::shared_ptr< OneDimensionalInterpolator< double, basic_mathematics::Vector6d > > interpolator =
            interpolators::createOneDimensionalInterpolator(
                timeHistoryOfState, interpolatorSettings );

    // Create ephemeris and return.
    return boost::make_shared< TabulatedCartesianEphemeris< > >(
                interpolator, observerName, referenceFrameName );
}
#endif

//! Function to create a ephemeris model.
boost::shared_ptr< ephemerides::Ephemeris > createBodyEphemeris(
        const boost::shared_ptr< EphemerisSettings > ephemerisSettings,
        const std::string& bodyName )
{
    // Declare return object.
    boost::shared_ptr< ephemerides::Ephemeris > ephemeris;

    // Check which type of ephemeris model is to be created.
    switch( ephemerisSettings->getEphemerisType( ) )
    {
    #if USE_CSPICE
    case direct_spice_ephemeris:
    {
        // Check consistency of type and class.
        boost::shared_ptr< DirectSpiceEphemerisSettings > directEphemerisSettings =
                boost::dynamic_pointer_cast< DirectSpiceEphemerisSettings >( ephemerisSettings );
        if( directEphemerisSettings == NULL )
        {
            throw std::runtime_error(
                        "Error, expected direct spice ephemeris settings for body " + bodyName );
        }
        else
        {
            // Create corresponding ephemeris object.
            ephemeris = boost::make_shared< SpiceEphemeris >(
                        bodyName,
                        directEphemerisSettings->getFrameOrigin( ),
                        directEphemerisSettings->getCorrectForStellarAbberation( ),
                        directEphemerisSettings->getCorrectForLightTimeAbberation( ),
                        directEphemerisSettings->getConvergeLighTimeAbberation( ),
                        directEphemerisSettings->getFrameOrientation( ) );
        }
        break;
    }
    case interpolated_spice:
    {
        // Check consistency of type and class.
        boost::shared_ptr< InterpolatedSpiceEphemerisSettings > interpolatedEphemerisSettings =
                boost::dynamic_pointer_cast< InterpolatedSpiceEphemerisSettings >(
                    ephemerisSettings );
        if( interpolatedEphemerisSettings == NULL )
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
                std::cerr<<"Warning, position of "<<bodyName<<" taken as barycenter of that body's "
                        <<"planetary system."<<std::endl;
            }

            // Create corresponding ephemeris object.
            ephemeris = createTabulatedEphemerisFromSpice(
                        inputName,
                        interpolatedEphemerisSettings->getInitialTime( ),
                        interpolatedEphemerisSettings->getFinalTime( ),
                        interpolatedEphemerisSettings->getTimeStep( ),
                        interpolatedEphemerisSettings->getFrameOrigin( ),
                        interpolatedEphemerisSettings->getFrameOrientation( ),
                        interpolatedEphemerisSettings->getInterpolatorSettings( ) );
        }
        break;
    }
    #endif
    case tabulated_ephemeris:
    {
        // Check consistency of type and class.
        boost::shared_ptr< TabulatedEphemerisSettings > tabulatedEphemerisSettings =
                boost::dynamic_pointer_cast< TabulatedEphemerisSettings >( ephemerisSettings );
        if( tabulatedEphemerisSettings == NULL )
        {
            throw std::runtime_error(
                        "Error, expected tabulated ephemeris settings for body " + bodyName );
        }
        else
        {
            // Create corresponding ephemeris object.
            ephemeris = boost::make_shared< TabulatedCartesianEphemeris< > >(
                        boost::make_shared<
                        interpolators::LagrangeInterpolator< double, basic_mathematics::Vector6d > >
                        ( tabulatedEphemerisSettings->getBodyStateHistory( ), 6,
                          interpolators::huntingAlgorithm,
                          interpolators::lagrange_cubic_spline_boundary_interpolation ),
                        tabulatedEphemerisSettings->getFrameOrigin( ),
                        tabulatedEphemerisSettings->getFrameOrientation( ) );
        }
        break;
    }
    case constant_ephemeris:
    {
        // Check consistency of type and class.
        boost::shared_ptr< ConstantEphemerisSettings > constantEphemerisSettings =
                boost::dynamic_pointer_cast< ConstantEphemerisSettings >( ephemerisSettings );
        if( constantEphemerisSettings == NULL )
        {
            throw std::runtime_error( "Error, expected constant ephemeris settings for " + bodyName );
        }
        else
        {
            // Create ephemeris
            ephemeris = boost::make_shared< ConstantEphemeris >(
                        boost::lambda::constant( constantEphemerisSettings->getConstantState( ) ),
                        constantEphemerisSettings->getFrameOrigin( ),
                        constantEphemerisSettings->getFrameOrientation( ) );
        }
        break;
    }
    case kepler_ephemeris:
    {
        // Check consistency of type and class.
        boost::shared_ptr< KeplerEphemerisSettings > keplerEphemerisSettings =
                boost::dynamic_pointer_cast< KeplerEphemerisSettings >( ephemerisSettings );
        if( keplerEphemerisSettings == NULL )
        {
            throw std::runtime_error( "Error, expected Kepler ephemeris settings for " + bodyName );
        }
        else
        {
            // Create ephemeris
            ephemeris = boost::make_shared< KeplerEphemeris >(
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
        boost::shared_ptr< ApproximatePlanetPositionSettings > approximateEphemerisSettings =
                boost::dynamic_pointer_cast< ApproximatePlanetPositionSettings >(
                    ephemerisSettings );
        if( approximateEphemerisSettings == NULL )
        {
            throw std::runtime_error(
                        "Error, expected approximate ephemeris settings for body " + bodyName );
        }
        else
        {
            // Create corresponding ephemeris object.
            if( approximateEphemerisSettings->getUseCircularCoplanarApproximation( ) )
            {
                ephemeris = boost::make_shared< ApproximatePlanetPositionsCircularCoplanar >(
                            approximateEphemerisSettings->getBodyIdentifier( ) );
            }
            else
            {
                ephemeris = boost::make_shared< ApproximatePlanetPositions >(
                            approximateEphemerisSettings->getBodyIdentifier( ) );
            }
        }
        break;
    }
    default:
    {
        throw std::runtime_error(
                    "Error, did not recognize ephemeris model settings type " +
                    boost::lexical_cast< std::string >( ephemerisSettings->getEphemerisType( ) ) );
    }
    }
    return ephemeris;

}

} // namespace simulation_setup

} // namespace tudat
