/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <vector>

#include "tudat/basics/testMacros.h"
#include "tudat/basics/utilities.h"
#include "tudat/io/matrixTextFileReader.h"

#include "tudat/math/interpolators/createInterpolator.h"
#include "tudat/io/basicInputOutput.h"

#include "tudat/interface/spice/spiceInterface.h"

namespace tudat
{
namespace unit_tests
{

using namespace interpolators;

BOOST_AUTO_TEST_SUITE( test_interpolator_vector_conversion )

BOOST_AUTO_TEST_CASE( testInterpolatorVectorConversion )
{
    spice_interface::loadStandardSpiceKernels( );

    std::map< double, Eigen::Vector6d > vector6dInterpolatorInput;
    std::map< double, Eigen::VectorXd > vectorXdInterpolatorInput;
    std::map< double, Eigen::MatrixXd > matrixXdInterpolatorInput;
    for( int i = 0; i < 30; i++ )
    {
        double time = static_cast< double >( i ) * 86400.0;
        vector6dInterpolatorInput[ time ] = spice_interface::getBodyCartesianStateAtEpoch(
                    "Moon", "Earth", "J2000", "None", time );
        vectorXdInterpolatorInput[ time ] = vector6dInterpolatorInput[ time ];
        matrixXdInterpolatorInput[ time ] = vector6dInterpolatorInput[ time ];
    }

    std::vector< double > interpolationTimes = { 1.0, 5.0 * 86400.0, 12.43 * 86400.0 };

    for( int i = 0; i < 4; i++ )
    {
        std::shared_ptr< InterpolatorSettings > interpolatorSettings;
        switch( i )
        {
        case 0:
            interpolatorSettings =
                    std::make_shared< InterpolatorSettings >( linear_interpolator );
            break;
        case 1:
            interpolatorSettings =
                    std::make_shared< InterpolatorSettings >( cubic_spline_interpolator );
            break;
        case 2:
            interpolatorSettings =
                    std::make_shared< LagrangeInterpolatorSettings >( 8 );
            break;
//        case 3:
//            interpolatorSettings =
//                    std::make_shared< InterpolatorSettings >( hermite_spline_interpolator );
//            break;
        case 3:
            interpolatorSettings =
                    std::make_shared< InterpolatorSettings >( piecewise_constant_interpolator );
            break;
        }
        std::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector6d > > direct6dInterpolator =
                createOneDimensionalInterpolator( vector6dInterpolatorInput, interpolatorSettings );
        std::shared_ptr< OneDimensionalInterpolator< double, Eigen::VectorXd > > directXdInterpolator =
                createOneDimensionalInterpolator( vectorXdInterpolatorInput, interpolatorSettings );
        std::shared_ptr< OneDimensionalInterpolator< double, Eigen::MatrixXd > > directXdMatrixInterpolator =
                createOneDimensionalInterpolator( matrixXdInterpolatorInput, interpolatorSettings );

        std::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector6d > > converted6dInterpolator =
                convertBetweenStaticDynamicEigenTypeInterpolators<  double, double, -1, 1, 6, 1 >( directXdInterpolator );
        std::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector6d > > converted6dInterpolator2 =
                convertBetweenStaticDynamicEigenTypeInterpolators<  double, double, -1, -1, 6, 1 >( directXdMatrixInterpolator );

        std::shared_ptr< OneDimensionalInterpolator< double, Eigen::VectorXd > > convertedXdInterpolator =
                convertBetweenStaticDynamicEigenTypeInterpolators<  double, double, 6, 1, -1, 1 >( direct6dInterpolator );
        std::shared_ptr< OneDimensionalInterpolator< double, Eigen::VectorXd > > convertedXdInterpolator2 =
                convertBetweenStaticDynamicEigenTypeInterpolators<  double, double, -1, -1, -1, 1 >( directXdMatrixInterpolator );

        std::shared_ptr< OneDimensionalInterpolator< double, Eigen::MatrixXd > > convertedXdMatrixInterpolator =
                convertBetweenStaticDynamicEigenTypeInterpolators<  double, double, 6, 1, -1, -1 >( direct6dInterpolator );
        std::shared_ptr< OneDimensionalInterpolator< double, Eigen::MatrixXd > > convertedXdMatrixInterpolator2 =
                convertBetweenStaticDynamicEigenTypeInterpolators<  double, double, -1, 1, -1, -1 >( directXdInterpolator );

        for( unsigned int j = 0; j < interpolationTimes.size( ); j++ )
        {
            std::vector< Eigen::MatrixXd > testMatrices;
            testMatrices.push_back( direct6dInterpolator->interpolate( interpolationTimes.at( j ) ) );
            testMatrices.push_back( directXdInterpolator->interpolate( interpolationTimes.at( j ) ) );
            testMatrices.push_back( directXdMatrixInterpolator->interpolate( interpolationTimes.at( j ) ) );

            testMatrices.push_back( converted6dInterpolator->interpolate( interpolationTimes.at( j ) ) );
            testMatrices.push_back( converted6dInterpolator2->interpolate( interpolationTimes.at( j ) ) );

            testMatrices.push_back( convertedXdInterpolator->interpolate( interpolationTimes.at( j ) ) );
            testMatrices.push_back( convertedXdInterpolator2->interpolate( interpolationTimes.at( j ) ) );

            testMatrices.push_back( convertedXdMatrixInterpolator->interpolate( interpolationTimes.at( j ) ) );
            testMatrices.push_back( convertedXdMatrixInterpolator2->interpolate( interpolationTimes.at( j ) ) );

            for( unsigned int k = 1; k < testMatrices.size( ); k++ )
            {
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            ( testMatrices.at( 0 ) ), ( testMatrices.at( k ) ), std::numeric_limits< double >::epsilon( ) );
            }

        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
