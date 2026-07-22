/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <memory>
#include <sstream>
#include <typeinfo>
#include <vector>

#include <boost/test/unit_test.hpp>
#include <Eigen/Core>

#include <cereal/archives/binary.hpp>

#include "tudat/astro/orbit_determination/podInputOutputTypes.h"
#include "tudat/io/serialization/registrations_estimation.h"

namespace tudat
{
namespace unit_tests
{

namespace
{

template< typename OutputType >
std::shared_ptr< OutputType > roundTripSerialize( const std::shared_ptr< OutputType >& original )
{
    std::stringstream serializationStream;

    {
        cereal::BinaryOutputArchive outputArchive( serializationStream );
        outputArchive( original );
    }

    std::shared_ptr< OutputType > roundTripped;

    {
        cereal::BinaryInputArchive inputArchive( serializationStream );
        inputArchive( roundTripped );
    }

    return roundTripped;
}

template< typename TimeType >
std::shared_ptr< simulation_setup::CovarianceAnalysisOutput< double, TimeType > > createCovarianceAnalysisOutput( )
{
    const Eigen::MatrixXd normalizedDesignMatrix( ( Eigen::MatrixXd( 3, 2 ) << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 ).finished( ) );
    const Eigen::VectorXd weightsMatrixDiagonal( ( Eigen::VectorXd( 3 ) << 7.0, 8.0, 9.0 ).finished( ) );
    const Eigen::VectorXd designMatrixTransformationDiagonal( ( Eigen::VectorXd( 2 ) << 10.0, 11.0 ).finished( ) );
    const Eigen::MatrixXd inverseNormalizedCovarianceMatrix( ( Eigen::MatrixXd( 2, 2 ) << 12.0, 13.0, 14.0, 15.0 ).finished( ) );
    const Eigen::MatrixXd normalizedDesignMatrixConsiderParameters( ( Eigen::MatrixXd( 3, 1 ) << 16.0, 17.0, 18.0 ).finished( ) );
    const Eigen::VectorXd considerNormalizationFactors( ( Eigen::VectorXd( 1 ) << 19.0 ).finished( ) );
    const Eigen::MatrixXd considerCovarianceContribution( ( Eigen::MatrixXd( 2, 2 ) << 20.0, 21.0, 22.0, 23.0 ).finished( ) );
    const Eigen::MatrixXd considerCovariance( ( Eigen::MatrixXd( 2, 2 ) << 24.0, 25.0, 26.0, 27.0 ).finished( ) );

    return std::make_shared< simulation_setup::CovarianceAnalysisOutput< double, TimeType > >( normalizedDesignMatrix,
                                                                                               weightsMatrixDiagonal,
                                                                                               designMatrixTransformationDiagonal,
                                                                                               inverseNormalizedCovarianceMatrix,
                                                                                               normalizedDesignMatrixConsiderParameters,
                                                                                               considerNormalizationFactors,
                                                                                               considerCovarianceContribution,
                                                                                               considerCovariance,
                                                                                               true );
}

template< typename TimeType >
std::shared_ptr< simulation_setup::EstimationOutput< double, TimeType > > createEstimationOutput( )
{
    const Eigen::Matrix< double, Eigen::Dynamic, 1 > parameterEstimate( ( Eigen::Vector2d( ) << 28.0, 29.0 ).finished( ) );
    const Eigen::VectorXd residuals( ( Eigen::Vector3d( ) << 30.0, 31.0, 32.0 ).finished( ) );
    const Eigen::MatrixXd normalizedDesignMatrix( ( Eigen::MatrixXd( 3, 2 ) << 33.0, 34.0, 35.0, 36.0, 37.0, 38.0 ).finished( ) );
    const Eigen::VectorXd weightsMatrixDiagonal( ( Eigen::VectorXd( 3 ) << 39.0, 40.0, 41.0 ).finished( ) );
    const Eigen::VectorXd designMatrixTransformationDiagonal( ( Eigen::VectorXd( 2 ) << 42.0, 43.0 ).finished( ) );
    const Eigen::MatrixXd inverseNormalizedCovarianceMatrix( ( Eigen::MatrixXd( 2, 2 ) << 44.0, 45.0, 46.0, 47.0 ).finished( ) );
    const double residualStandardDeviation = 48.0;
    const int bestIteration = 1;
    const std::vector< Eigen::VectorXd > residualHistory = { ( Eigen::Vector3d( ) << 49.0, 50.0, 51.0 ).finished( ),
                                                             ( Eigen::Vector3d( ) << 52.0, 53.0, 54.0 ).finished( ) };
    const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > parameterHistory = { ( Eigen::Vector2d( ) << 55.0, 56.0 ).finished( ),
                                                                                         ( Eigen::Vector2d( ) << 57.0, 58.0 ).finished( ) };
    const Eigen::MatrixXd normalizedDesignMatrixConsiderParameters( ( Eigen::MatrixXd( 3, 1 ) << 59.0, 60.0, 61.0 ).finished( ) );
    const Eigen::VectorXd considerNormalizationFactors( ( Eigen::VectorXd( 1 ) << 62.0 ).finished( ) );
    const Eigen::MatrixXd covarianceConsiderContribution( ( Eigen::MatrixXd( 2, 2 ) << 63.0, 64.0, 65.0, 66.0 ).finished( ) );
    const Eigen::MatrixXd considerCovariance( ( Eigen::MatrixXd( 2, 2 ) << 67.0, 68.0, 69.0, 70.0 ).finished( ) );

    return std::make_shared< simulation_setup::EstimationOutput< double, TimeType > >( parameterEstimate,
                                                                                       residuals,
                                                                                       normalizedDesignMatrix,
                                                                                       weightsMatrixDiagonal,
                                                                                       designMatrixTransformationDiagonal,
                                                                                       inverseNormalizedCovarianceMatrix,
                                                                                       residualStandardDeviation,
                                                                                       bestIteration,
                                                                                       residualHistory,
                                                                                       parameterHistory,
                                                                                       normalizedDesignMatrixConsiderParameters,
                                                                                       considerNormalizationFactors,
                                                                                       covarianceConsiderContribution,
                                                                                       considerCovariance,
                                                                                       true,
                                                                                       false );
}

}  // namespace

BOOST_AUTO_TEST_SUITE( test_PodInputOutputTypes_serialization )

BOOST_AUTO_TEST_CASE( test_CovarianceAnalysisOutputSerialization )
{
    std::vector< std::shared_ptr< simulation_setup::CovarianceAnalysisOutput< double, double > > > resultsVector = {
        createCovarianceAnalysisOutput< double >( ), createEstimationOutput< double >( )
    };

    for( const auto& result : resultsVector )
    {
        auto roundTripped = roundTripSerialize( result );
        BOOST_REQUIRE( roundTripped != nullptr );
        BOOST_CHECK_MESSAGE( *result == *roundTripped,
                             "Round-tripped object is not equal to original. Type: " << typeid( *result ).name( ) );
    }
    // checkCovarianceAnalysisOutputRoundTrip< double >( );
    // checkCovarianceAnalysisOutputRoundTrip< tudat::Time >( );

    // checkEstimationOutputRoundTrip< double >( );
    // checkEstimationOutputRoundTrip< tudat::Time >( );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
