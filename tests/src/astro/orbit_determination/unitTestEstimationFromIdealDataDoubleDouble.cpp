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


#include <limits>

#include <boost/test/unit_test.hpp>

#include "tudat/simulation/estimation_setup/orbitDeterminationTestCases.h"


namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_estimation_from_positions )



//! This test checks, for double states/observables and double time, if the orbit determination correctly converges
//! when simulating data, perturbing the dynamical parameters, and then retrieving the original parameters
BOOST_AUTO_TEST_CASE( test_EstimationFromPosition )
{
    for( int simulationType = 0; simulationType < 6; simulationType++ )
    {

        std::cout << "=============================================== Running Case: " << simulationType << std::endl;

        // Simulate estimated parameter error.
        Eigen::VectorXd totalError;

        totalError = executePlanetaryParameterEstimation< double, double >( simulationType ).second;

        // Adjust tolerance based on simulation settings
        double toleranceMultiplier = 20.0;

        // Check error.
        for( unsigned int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( totalError( j ), toleranceMultiplier * 5.0E-3 );
        }

        for( unsigned int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( totalError( j + 3 ), toleranceMultiplier * 1.0E-7 );
        }

        BOOST_CHECK_SMALL( totalError( 6 ), toleranceMultiplier * 1.0E3 );
        std::cout <<"Total error: "<< totalError.transpose( ) << std::endl;
    }

    std::pair< std::shared_ptr< simulation_setup::PodOutput< double > >,
    std::shared_ptr< simulation_setup::PodInput< double, double > > > podDataOutput;
    Eigen::VectorXd estimationError = tudat::unit_tests::executeEarthOrbiterParameterEstimation< double, double >(
                 podDataOutput );

    for( unsigned int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( estimationError( i ) ), 1.0E-5 );
        BOOST_CHECK_SMALL( std::fabs( estimationError( i + 3 ) ), 1.0E-8 );
        BOOST_CHECK_SMALL( std::fabs( estimationError( i + 18 ) ), 1.0E-4 );
    }

    BOOST_CHECK_SMALL( std::fabs( estimationError( 6 ) ), 5.0E-7 );
    BOOST_CHECK_SMALL( std::fabs( estimationError( 7 ) ), 5.0E-7 );
    BOOST_CHECK_SMALL( std::fabs( estimationError( 8 ) ), 5.0E-7 );

    BOOST_CHECK_SMALL( std::fabs( estimationError( 9 ) ), 5.0E-13 );
    BOOST_CHECK_SMALL( std::fabs( estimationError( 10 ) ), 5.0E-13 );
    BOOST_CHECK_SMALL( std::fabs( estimationError( 11 ) ), 5.0E-13 );

    for( unsigned int i = 0; i < 6; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( estimationError( i + 12 ) ), 5.0E-13 );
    }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}


