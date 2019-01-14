/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <Tudat/SimulationSetup/tudatEstimationHeader.h>
#include "Tudat/SimulationSetup/PropagationSetup/fullPropagationLambertTargeter.h"
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <Eigen/Core>
#include "Tudat/Basics/testMacros.h"


namespace tudat
{

namespace unit_tests
{

using namespace tudat;

//! Test if the difference between the Lambert targeter solution and the full dynamics problem is computed correctly.
BOOST_AUTO_TEST_CASE( testFullPropagationLambertTargeter )
{

    std::cout.precision(20);

    double initialTime = 0.0;
    double fixedStepSize = 1.0;

    Eigen::Vector3d cartesianPositionAtDeparture ( 2.0 * 6.378136e6, 0.0, 0.0 );
    Eigen::Vector3d cartesianPositionAtArrival ( 2.0 * 6.378136e6, 2.0 * std::sqrt( 3.0 ) * 6.378136e6, 0.0 );

    double timeOfFlight = 806.78 * 5.0;

    std::vector< std::string > bodiesToPropagate; bodiesToPropagate.push_back("spacecraft");
    std::vector< std::string > centralBodies; centralBodies.push_back("Earth");

    // Define integrator settings.
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared < numerical_integrators::IntegratorSettings < > >
                ( numerical_integrators::rungeKutta4, initialTime, fixedStepSize);

    // Define the body map.
    simulation_setup::NamedBodyMap bodyMap = propagators::setupBodyMapLambertTargeter("Earth", "spacecraft");
    basic_astrodynamics::AccelerationMap accelerationModelMap = propagators::setupAccelerationMapLambertTargeter(
                "Earth", "spacecraft", bodyMap);

    std::vector< std::string > departureAndArrivalBodies;
    departureAndArrivalBodies.push_back("Earth");
    departureAndArrivalBodies.push_back("Mars");

   // Compute the difference in state between the full problem and the Lambert targeter solution at departure and at arrival
    std::pair< Eigen::Vector6d, Eigen::Vector6d > differenceState =
            propagators::getDifferenceFullPropagationWrtLambertTargeterAtDepartureAndArrival(cartesianPositionAtDeparture,
             cartesianPositionAtArrival, timeOfFlight, bodyMap, accelerationModelMap, bodiesToPropagate,
             centralBodies, integratorSettings, departureAndArrivalBodies);

    Eigen::Vector6d differenceStateAtDeparture = differenceState.first;
    Eigen::Vector6d differenceStateAtArrival = differenceState.second;

    std::cout << "differenceStateAtDeparture: " << differenceStateAtDeparture << "\n\n";
    std::cout << "differenceStateAtArrival: " << differenceStateAtArrival << "\n\n";


    for( int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( differenceStateAtDeparture( i ) ), 1.0 );
        BOOST_CHECK_SMALL( std::fabs( differenceStateAtDeparture( i + 3 ) ), 1.0E-6 );
        BOOST_CHECK_SMALL( std::fabs( differenceStateAtArrival( i ) ), 1.0 );
        BOOST_CHECK_SMALL( std::fabs( differenceStateAtArrival( i + 3 ) ), 1.0E-6 );
    }

}

}

}


