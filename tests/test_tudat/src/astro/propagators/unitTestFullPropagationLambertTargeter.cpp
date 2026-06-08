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

BOOST_AUTO_TEST_SUITE( testFullPropagationLambertTargeter )

//! Test if the difference between the Lambert targeter solution and the full dynamics problem is computed correctly.
BOOST_AUTO_TEST_CASE( testFullPropagationLambertTargeterBasic )
{
    std::cout.precision( 20 );

    double initialTime = 0.0;

    Eigen::Vector3d cartesianPositionAtDeparture( 2.0 * 6.378136e6, 0.0, 0.0 );
    Eigen::Vector3d cartesianPositionAtArrival( 2.0 * 6.378136e6, 2.0 * std::sqrt( 3.0 ) * 6.378136e6, 0.0 );

    double timeOfFlight = 806.78 * 5.0;
    double fixedStepSize = timeOfFlight / 10000.0;

    std::string bodyToPropagate = "spacecraft";
    std::string centralBody = "Earth";

    // Define integrator settings.
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings<> >(
                    numerical_integrators::rungeKutta4, initialTime, fixedStepSize );

    std::vector< std::string > departureAndArrivalBodyNames = { "departure", "arrival" };
    std::pair< std::string, std::string > departureAndArrivalBodies =
            std::make_pair( departureAndArrivalBodyNames.at( 0 ), departureAndArrivalBodyNames.at( 1 ) );

    // Define the system of bodies.
    simulation_setup::SystemOfBodies bodies = propagators::setupBodyMapFromUserDefinedStatesForLambertTargeter(
            "Earth", "spacecraft", departureAndArrivalBodyNames, cartesianPositionAtDeparture, cartesianPositionAtArrival );

    basic_astrodynamics::AccelerationMap accelerationModelMap =
            propagators::setupAccelerationMapLambertTargeter( "Earth", "spacecraft", bodies );

    std::map< double, Eigen::Vector6d > lambertTargeterResult;
    std::map< double, Eigen::Vector6d > fullProblemResult;
    std::map< double, Eigen::VectorXd > dependentVariableResult;
    propagateLambertTargeterAndFullProblem( timeOfFlight,
                                            initialTime,
                                            bodies,
                                            accelerationModelMap,
                                            bodyToPropagate,
                                            centralBody,
                                            departureAndArrivalBodies,
                                            integratorSettings,
                                            lambertTargeterResult,
                                            fullProblemResult,
                                            dependentVariableResult,
                                            false );

    for( auto stateIterator : lambertTargeterResult )
    {
        for( int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL(
                    std::fabs( lambertTargeterResult.at( stateIterator.first )( i ) - fullProblemResult.at( stateIterator.first )( i ) ),
                    1.0E-5 );
            BOOST_CHECK_SMALL( std::fabs( lambertTargeterResult.at( stateIterator.first )( i + 3 ) -
                                          fullProblemResult.at( stateIterator.first )( i + 3 ) ),
                               1.0E-9 );
        }
    }
}

}  // namespace unit_tests

}  // namespace tudat
