/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>

#include "tudat/io/trackingData.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_tracking_data )

data::TrackingData< double, double > createVectorTrackingData( )
{
    std::vector< Eigen::VectorXd > observations( 2 );
    observations[ 0 ] = Eigen::Vector2d( 1.0, 2.0 );
    observations[ 1 ] = Eigen::Vector2d( 3.0, 4.0 );

    data::PlainLinkDefinition linkEnds;
    const std::vector< double > epochs = { 10.0, 20.0 };
    return data::TrackingData< double, double >( "angular_position", linkEnds, observations, epochs, "receiver" );
}

BOOST_AUTO_TEST_CASE( testConcatenatedObservationsForVectorObservables )
{
    data::TrackingData< double, double > trackingData = createVectorTrackingData( );

    const Eigen::VectorXd concatenatedObservations = trackingData.getObservationsVector( );

    BOOST_REQUIRE_EQUAL( concatenatedObservations.rows( ), 4 );
    BOOST_CHECK_EQUAL( concatenatedObservations( 0 ), 1.0 );
    BOOST_CHECK_EQUAL( concatenatedObservations( 1 ), 2.0 );
    BOOST_CHECK_EQUAL( concatenatedObservations( 2 ), 3.0 );
    BOOST_CHECK_EQUAL( concatenatedObservations( 3 ), 4.0 );
}

BOOST_AUTO_TEST_CASE( testConcatenatedObservationEpochsForVectorObservables )
{
    data::TrackingData< double, double > trackingData = createVectorTrackingData( );

    // One epoch is stored per vector-valued observation. The flattened form repeats
    // that epoch for each observable component so it aligns with the observation vector.
    const Eigen::VectorXd concatenatedEpochs = trackingData.getObservationEpochsVector( );

    BOOST_REQUIRE_EQUAL( concatenatedEpochs.rows( ), 4 );
    BOOST_CHECK_EQUAL( concatenatedEpochs( 0 ), 10.0 );
    BOOST_CHECK_EQUAL( concatenatedEpochs( 1 ), 10.0 );
    BOOST_CHECK_EQUAL( concatenatedEpochs( 2 ), 20.0 );
    BOOST_CHECK_EQUAL( concatenatedEpochs( 3 ), 20.0 );
}

BOOST_AUTO_TEST_CASE( testConcatenatedWeightsForVectorObservables )
{
    data::TrackingData< double, double > trackingData = createVectorTrackingData( );

    std::vector< Eigen::VectorXd > weights( 2 );
    weights[ 0 ] = Eigen::Vector2d( 0.1, 0.2 );
    weights[ 1 ] = Eigen::Vector2d( 0.3, 0.4 );
    trackingData.setObservationWeights( weights );

    const Eigen::VectorXd concatenatedWeights = trackingData.getObservationWeightsVector( );

    BOOST_REQUIRE_EQUAL( concatenatedWeights.rows( ), 4 );
    BOOST_CHECK_EQUAL( concatenatedWeights( 0 ), 0.1 );
    BOOST_CHECK_EQUAL( concatenatedWeights( 1 ), 0.2 );
    BOOST_CHECK_EQUAL( concatenatedWeights( 2 ), 0.3 );
    BOOST_CHECK_EQUAL( concatenatedWeights( 3 ), 0.4 );
}

BOOST_AUTO_TEST_CASE( testVectorObservableSizeValidation )
{
    data::TrackingData< double, double > trackingData = createVectorTrackingData( );

    std::vector< Eigen::VectorXd > wrongNumberOfWeights( 1 );
    wrongNumberOfWeights[ 0 ] = Eigen::Vector2d( 0.1, 0.2 );
    BOOST_CHECK_THROW( trackingData.setObservationWeights( wrongNumberOfWeights ), std::runtime_error );

    std::vector< Eigen::VectorXd > wrongSingleWeightSize( 2 );
    wrongSingleWeightSize[ 0 ] = Eigen::Vector3d( 0.1, 0.2, 0.3 );
    wrongSingleWeightSize[ 1 ] = Eigen::Vector3d( 0.4, 0.5, 0.6 );
    BOOST_CHECK_THROW( trackingData.setObservationWeights( wrongSingleWeightSize ), std::runtime_error );

    std::vector< Eigen::VectorXd > wrongNumberOfCorrections( 1 );
    wrongNumberOfCorrections[ 0 ] = Eigen::Vector2d( 0.1, 0.2 );
    BOOST_CHECK_THROW( trackingData.setObservationCorrections( wrongNumberOfCorrections ), std::runtime_error );

    std::vector< Eigen::VectorXd > wrongSingleCorrectionSize( 2 );
    wrongSingleCorrectionSize[ 0 ] = Eigen::Vector3d( 0.1, 0.2, 0.3 );
    wrongSingleCorrectionSize[ 1 ] = Eigen::Vector3d( 0.4, 0.5, 0.6 );
    BOOST_CHECK_THROW( trackingData.setObservationCorrections( wrongSingleCorrectionSize ), std::runtime_error );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
