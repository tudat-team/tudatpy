/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <tudat/simulation/estimation.h>
#include "tudat/simulation/propagation_setup/propagationLambertTargeterFullProblem.h"
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <Eigen/Core>
#include "tudat/basics/testMacros.h"

#include "tudat/astro/ephemerides/approximatePlanetPositions.h"
#include "tudat/astro/trajectory_design/trajectory.h"
#include "tudat/astro/trajectory_design/exportTrajectory.h"
#include "tudat/astro/trajectory_design/planetTrajectory.h"

using namespace tudat;
using namespace tudat::propagators;

namespace tudat
{

namespace unit_tests
{



//! Function to setup a body map corresponding to the assumptions of the Lambert targeter,
//! using default ephemerides for the central body only, while the positions of departure and arrival bodies are provided as inputs.
/*!
 * Function to setup a Lambert targeter map. The body map contains the central, departure and arrival bodies and the body to be propagated.
 * The positions of the departure and arrival bodies are defined by the user and provided as inputs.
 * \param nameCentralBody Name of the central body.
 * \param nameBodyToPropagate Name of the body to be propagated.
 * \param departureAndArrivalBodies Vector containing the names of the departure and arrival bodies.
 * \param cartesianPositionAtDeparture Vector containing the position coordinates of the departure body [m].
 * \param cartesianPositionAtArrival Vector containing the position coordinates of the arrival body [m].
 * \return Body map for the Lambert targeter.
 */
simulation_setup::NamedBodyMap setupBodyMapFromUserDefinedStatesForLambertTargeter(
        const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate,
        const std::pair< std::string, std::string >&  departureAndArrivalBodies,
        const Eigen::Vector3d& cartesianPositionAtDeparture,
        const Eigen::Vector3d& cartesianPositionAtArrival )
{
    spice_interface::loadStandardSpiceKernels( );

    std::string frameOrigin = "SSB";
    std::string frameOrientation = "ECLIPJ2000";

    // Create ephemeris vector for departure and arrival bodies.
    std::vector< ephemerides::EphemerisPointer > ephemerisVectorDepartureAndArrivalBodies(2);
    ephemerisVectorDepartureAndArrivalBodies[ 0 ] = std::make_shared< ephemerides::ConstantEphemeris > (
                ( Eigen::Vector6d( ) << cartesianPositionAtDeparture,
                  0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );
    ephemerisVectorDepartureAndArrivalBodies[ 1 ] = std::make_shared< ephemerides::ConstantEphemeris >(
                ( Eigen::Vector6d( ) << cartesianPositionAtArrival,
                  0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );

    // Create body map.
    simulation_setup::NamedBodyMap bodyMap =
            propagators::setupBodyMapFromUserDefinedEphemeridesForLambertTargeter(
                nameCentralBody, nameBodyToPropagate, departureAndArrivalBodies, ephemerisVectorDepartureAndArrivalBodies );

    return bodyMap;

}

BOOST_AUTO_TEST_SUITE( testFullPropagationLambertTargeter )

//! Test if the difference between the Lambert targeter solution and the full dynamics problem is computed correctly.
BOOST_AUTO_TEST_CASE( testFullPropagationLambertTargeterBasic )
{

    std::cout.precision(20);

    double initialTime = 0.0;

    Eigen::Vector3d cartesianPositionAtDeparture ( 2.0 * 6.378136e6, 0.0, 0.0 );
    Eigen::Vector3d cartesianPositionAtArrival ( 2.0 * 6.378136e6, 2.0 * std::sqrt( 3.0 ) * 6.378136e6, 0.0 );

    double timeOfFlight = 806.78 * 5.0;
    double fixedStepSize = timeOfFlight / 10000.0;

    std::string bodyToPropagate = "spacecraft" ;
    std::string centralBody = "Earth";

    // Define integrator settings.
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared < numerical_integrators::IntegratorSettings < > >
            ( numerical_integrators::rungeKutta4, initialTime, fixedStepSize);


    std::pair< std::string, std::string > departureAndArrivalBodies =
            std::make_pair( "departure", "arrival" );

<<<<<<< HEAD
    // Define the system of bodies.
    simulation_setup::SystemOfBodies bodies = propagators::setupBodyMapFromUserDefinedStatesForLambertTargeter("Earth", "spacecraft", departureAndArrivalBodies,
                                                                                      cartesianPositionAtDeparture, cartesianPositionAtArrival);

    basic_astrodynamics::AccelerationMap accelerationModelMap = propagators::setupAccelerationMapLambertTargeter(
                "Earth", "spacecraft", bodies);


   // Compute the difference in state between the full problem and the Lambert targeter solution at departure and at arrival
   std::pair< Eigen::Vector6d, Eigen::Vector6d > differenceState =
            propagators::getDifferenceFullPropagationWrtLambertTargeterAtDepartureAndArrival(cartesianPositionAtDeparture,
             cartesianPositionAtArrival, timeOfFlight, initialTime, bodies, accelerationModelMap, bodyToPropagate,
             centralBody, integratorSettings, departureAndArrivalBodies, false);
=======
    // Define the body map.
    simulation_setup::NamedBodyMap bodyMap = setupBodyMapFromUserDefinedStatesForLambertTargeter(
                "Earth", "spacecraft", departureAndArrivalBodies, cartesianPositionAtDeparture, cartesianPositionAtArrival );

    basic_astrodynamics::AccelerationMap accelerationModelMap = propagators::setupAccelerationMapLambertTargeter(
                "Earth", "spacecraft", bodyMap );
>>>>>>> dominic-origin/features/mission_segments_refactor

    std::map< double, Eigen::Vector6d > lambertTargeterResult;
    std::map< double, Eigen::Vector6d > fullProblemResult;
    std::map< double, Eigen::VectorXd > dependentVariableResult;
    propagateLambertTargeterAndFullProblem(
                timeOfFlight, initialTime, bodyMap, accelerationModelMap, bodyToPropagate,
                centralBody, departureAndArrivalBodies, integratorSettings, lambertTargeterResult, fullProblemResult,
                dependentVariableResult, false );

    for( auto stateIterator : lambertTargeterResult )
    {
        for( int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( std::fabs( lambertTargeterResult.at( stateIterator.first )( i ) -
                                          fullProblemResult.at( stateIterator.first )( i ) ), 1.0E-5 );
            BOOST_CHECK_SMALL( std::fabs( lambertTargeterResult.at( stateIterator.first )( i + 3 ) -
                                          fullProblemResult.at( stateIterator.first )( i + 3 ) ), 1.0E-9 );
        }
    }
}

}

}


