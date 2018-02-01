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


#include <limits>

#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/OrbitDetermination/UnitTests/orbitDeterminationTestCases.h"


namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_estimation_from_positions )



//! This test checks whether observation biases are correctly estimated, using a variety of different settings
//! for types of observables/biases.
BOOST_AUTO_TEST_CASE( test_EstimationFromPosition )
{
    for( int estimateRangeBiases = 0; estimateRangeBiases < 2; estimateRangeBiases++ )
    {
        for( int estimateTwoWayBiases = 0; estimateTwoWayBiases < 2; estimateTwoWayBiases++ )
        {
            for( int useSingleBiasModel = 0; useSingleBiasModel < 2; useSingleBiasModel++ )
            {
                for( int estimateAbsoluteBiases = 0; estimateAbsoluteBiases < 2; estimateAbsoluteBiases++ )
                {
                    for( int estimateMultiArcBiases = 0; estimateMultiArcBiases < 2; estimateMultiArcBiases++ )
                    {
                        std::cout << "=========== Running Case: " << estimateRangeBiases <<" "<<estimateTwoWayBiases<<" "<<
                                     useSingleBiasModel<<" "<<estimateAbsoluteBiases<<" "<<estimateMultiArcBiases<<std::endl;

                        // Simulate estimated parameter error.
                        Eigen::VectorXd totalError = executeEarthOrbiterBiasEstimation< double, double >(
                                    estimateRangeBiases, estimateTwoWayBiases, useSingleBiasModel, estimateAbsoluteBiases, false,
                                    estimateMultiArcBiases ).first;

                        for( unsigned int j = 0; j < 3; j++ )
                        {
                            BOOST_CHECK_SMALL( std::fabs( totalError( j ) ), 1.0E-5 );
                            BOOST_CHECK_SMALL( std::fabs( totalError( j + 3 ) ), 1.0E-8 );
                        }

                        for( unsigned int j = 6; j < totalError.rows( ); j++ )
                        {
                            if( estimateAbsoluteBiases )
                            {
                                if( estimateRangeBiases )
                                {
                                    BOOST_CHECK_SMALL( std::fabs( totalError( j ) ), 1.0E-7);
                                }
                                else
                                {
                                    BOOST_CHECK_SMALL( std::fabs( totalError( j ) ), 1.0E-18 );
                                }
                            }
                            else if( !estimateMultiArcBiases )
                            {
                                BOOST_CHECK_SMALL( std::fabs( totalError( j ) ), 1.0E-14 );
                            }
                            else
                            {
                                BOOST_CHECK_SMALL( std::fabs( totalError( j ) ), 1.0E-13 );
                            }
                        }
                    }

                }
            }
        }
    }

    BOOST_CHECK_EQUAL( executeEarthOrbiterBiasEstimation( true, false, true, true, true ).second, true );
}

BOOST_AUTO_TEST_SUITE_END( )

}

}


