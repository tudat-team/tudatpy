/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *        Wakker, K.F. "astro I, AE4-874", Delft University of Technology, 2007.
 *        Howell, K.C. Three-dimensional, periodic, 'Halo' orbits, Celestial Mechanics, 32. 53-71,
 *          1984.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/simulation/simulation.h"
#include "tudat/basics/testMacros.h"
#include "tudat/simulation/propagation_setup/propagationCR3BPFullProblem.h"


namespace tudat
{

namespace unit_tests
{

using namespace tudat;
using namespace propagators;
using namespace orbital_element_conversions;

//! Test if state derivative for circular restricted three-body problem is computed correctly.
BOOST_AUTO_TEST_CASE( testFullPropagationCircularRestrictedThreeBodyProblem )
{
    double initialTime = 0.0;
    double finalTime = 120000000.0;

    std::vector < std::string > bodiesCR3BP;
    bodiesCR3BP.push_back( "Sun" );
    bodiesCR3BP.push_back( "Earth" );

    simulation_setup::SystemOfBodies bodies = setupBodyMapCR3BP(
                physical_constants::ASTRONOMICAL_UNIT, "Sun", "Earth", "Spacecraft" );

    // Spacecraft properties
    bodies.at( "Spacecraft" )->setConstantBodyMass( 100.0 );

    // Initialization of the spacecraft state
    Eigen::Vector6d initialState;
    initialState[0] = 2.991957413820000e+10;
    initialState[1] = 1.295555563704656e+11;
    initialState[2] = 0.0;
    initialState[3] = -2.579433850734350e+04;
    initialState[4] = 5.956947312313238e+03;
    initialState[5] = 0.0;


    // Define propagator settings variables.
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;
    bodiesToPropagate.push_back( "Spacecraft" );
    centralBodies.push_back( "SSB" );

    basic_astrodynamics::AccelerationMap accelerationModelMap = setupAccelerationMapCR3BP(
                "Sun", "Earth", bodiesToPropagate.at( 0 ), centralBodies.at( 0 ), bodies );


    // Create integrator settings
    const double fixedStepSize = 1000.0;
    std::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings =
            std::make_shared < numerical_integrators::IntegratorSettings < > >
            ( numerical_integrators::rungeKutta4, initialTime, fixedStepSize );


    // calculate the difference between CR3BP and full problem
    Eigen::Vector6d stateDifference = getFinalStateDifferenceFullPropagationWrtCR3BP(
                initialTime, finalTime, initialState, integratorSettings, accelerationModelMap,
                bodiesToPropagate, centralBodies,bodies, bodiesCR3BP );

    for( int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( stateDifference( i ) ), 1.0 );
        BOOST_CHECK_SMALL( std::fabs( stateDifference( i + 3 ) ), 1.0E-6 );
    }

}

}

}
