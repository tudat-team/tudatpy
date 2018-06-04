/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
namespace tudat
{
namespace unit_tests
{

using namespace mathematical_constants;

BOOST_AUTO_TEST_SUITE( test_cr3bp_propagation )

//! Test if CR3BP propagation is working correctly
BOOST_AUTO_TEST_CASE( testCR3BPPropagation )
{
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::basic_mathematics;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;

    // Test three different cases of CR3BP propagation. Reference values independently generated bvy Y. Liu
    for( int testCase = 0; testCase < 3; testCase++ )
    {
        // Define initial normalized state for each case
        Eigen::Vector6d massLessBodyInitialCartesianState;
        if( testCase == 0 )
        {
            massLessBodyInitialCartesianState[ xCartesianPositionIndex ] = 0.994;
            massLessBodyInitialCartesianState[ yCartesianPositionIndex ] = 0.853;
            massLessBodyInitialCartesianState[ zCartesianPositionIndex ] = 0.312;
            massLessBodyInitialCartesianState[ xCartesianVelocityIndex ] = 0.195;
            massLessBodyInitialCartesianState[ yCartesianVelocityIndex ] = -0.211;
            massLessBodyInitialCartesianState[ zCartesianVelocityIndex ] = 0.15;
        }
        else if( testCase == 1 )
        {
            massLessBodyInitialCartesianState[ xCartesianPositionIndex ] = 1.05;
            massLessBodyInitialCartesianState[ yCartesianPositionIndex ] = 0.13;
            massLessBodyInitialCartesianState[ zCartesianPositionIndex ] = 0.0;
            massLessBodyInitialCartesianState[ xCartesianVelocityIndex ] = 0.2;
            massLessBodyInitialCartesianState[ yCartesianVelocityIndex ] = 0.2;
            massLessBodyInitialCartesianState[ zCartesianVelocityIndex ] = 0.0;
        }
        else if( testCase == 2 )
        {
            massLessBodyInitialCartesianState[ xCartesianPositionIndex ] = 0.988;
            massLessBodyInitialCartesianState[ yCartesianPositionIndex ] = 0.5;
            massLessBodyInitialCartesianState[ zCartesianPositionIndex ] = -0.1;
            massLessBodyInitialCartesianState[ xCartesianVelocityIndex ] = 0.4;
            massLessBodyInitialCartesianState[ yCartesianVelocityIndex ] = 0.1;
            massLessBodyInitialCartesianState[ zCartesianVelocityIndex ] = -0.1;
        }

        // Set normalized-mass parameter and propagation intervals
        double massParameter = 0.0;
        double simulationStartEpoch = 0.0;
        double simulationEndEpoch = 0.0;
        double timeStep = 0.0;
        if( testCase == 0 )
        {
            massParameter = 2.528e-5;
            simulationEndEpoch = 20.0;
            timeStep = 0.0001;
        }
        else if( testCase == 1 )
        {
            massParameter = 1.21506683e-2;
            simulationEndEpoch = 10.0;
            timeStep = 0.0001;
        }
        else if( testCase == 2 )
        {
            massParameter = 9.537e-4;
            simulationEndEpoch = 15.0;
            timeStep = 0.001;
        }

        // Set integrator settings
        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< IntegratorSettings< > >
                ( rungeKutta4, simulationStartEpoch, timeStep );

        // Propagate dynamics
        std::map< double, Eigen::Vector6d > stateHistory = performCR3BPIntegration(
                    integratorSettings, massParameter,
                    massLessBodyInitialCartesianState,
                    simulationEndEpoch );

        // Retrieve dynamics and set iterator at final time
        std::map< double, Eigen::Vector6d >::reverse_iterator stateIterator =
                stateHistory.rbegin( );

        // Check for overshoot, and increment iterator if overshoot found.
        if( stateIterator->first - simulationEndEpoch > timeStep / 2.0 )
        {
            stateIterator++;
        }

        // Set stets computed independently by Y. Loi
        Eigen::Vector6d expectedFinalState;
        if( testCase == 0 )
        {
            expectedFinalState << -1.34313636385140, -1.54200249942130, -0.416194453794142,  -0.863033291171519,
                    1.12530842202949, 0.181821699265344;

        }
        else if( testCase == 1 )
        {
            expectedFinalState << -0.636895857906707, 0.977767324698871, 0, 0.133701380037512, 0.0704414501083985, 0;

        }
        else if( testCase == 2 )
        {
            expectedFinalState << 0.770391451297307, -1.53581492635515, -0.181256821037128,  -0.783687153930728,
                    -0.805869201715479, -0.0200378316314463;

        }

        // Check Tudat output with reference values. Remaining differences can be attributed to different integrators
        for( unsigned int i = 0; i < 6; i++ )
        {
            BOOST_CHECK_SMALL( std::fabs( stateIterator->second( i ) - expectedFinalState( i ) ), 1.0E-10 );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}
