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


#include <limits>

#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/OrbitDetermination/UnitTests/orbitDeterminationTestCases.h"


namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_estimation_from_positions )



//! This test checks, for a variety of data types/floating point types, if the orbit determination correctly converges
//! when simulating data, perturbing the dynamical parameters, and then retrieving the original parameters
BOOST_AUTO_TEST_CASE( test_EstimationFromPosition )
{
    for( int simulationType = 0; simulationType < 4; simulationType++ )
    {
        for( unsigned int i = 0; i < 4; i++ )
        {
            std::cout<<"=============================================== Running Case: "<<i<<" "<<simulationType<<std::endl;

            // Simulate estimated parameter error.
            Eigen::VectorXd totalError;
            if( i == 0 )
            {
                totalError = executeParameterEstimation< double, double >( simulationType ).second;
            }
            else if( i == 1 )
            {
                totalError = executeParameterEstimation< double, long double >( simulationType ).second;
            }
            else if( i == 2 )
            {
                totalError = executeParameterEstimation< Time, double >( simulationType ).second;
            }
            else if( i == 3 )
            {
                totalError = executeParameterEstimation< Time, long double >( simulationType ).second;
            }

            // Adjust tolerance based on simulation settings
            double toleranceMultiplier = 1.0;
            if( i % 2 == 1 )
            {
                toleranceMultiplier *= 1.0E-3;

                if( simulationType > 0 )
                {
                    toleranceMultiplier *= 100.0;
                }
            }
            else if( simulationType > 0 )
            {
                toleranceMultiplier *= 20.0;
            }

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
            std::cout<<totalError.transpose( )<<std::endl;
        }
    }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}


