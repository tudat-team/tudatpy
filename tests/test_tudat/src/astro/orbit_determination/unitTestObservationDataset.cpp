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

#include <boost/test/unit_test.hpp>

#include <cmath>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation_setup/observationDataset.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::observation_models;

BOOST_AUTO_TEST_SUITE( test_observation_dataset )

LinkDefinition createOneWayLinkDefinition( const std::string& stationName )
{
    LinkEnds linkEnds;
    linkEnds[ transmitter ] = LinkEndId( "Earth", stationName );
    linkEnds[ receiver ] = LinkEndId( "Vehicle", "" );
    return LinkDefinition( linkEnds );
}

std::shared_ptr< SingleObservationSet< double, double > > createManualObservationSet( const ObservableType observableType,
                                                                                      const LinkDefinition& linkDefinition,
                                                                                      const std::vector< Eigen::VectorXd >& observations,
                                                                                      const std::vector< double >& times,
                                                                                      const std::vector< Eigen::VectorXd >& weights,
                                                                                      const std::vector< Eigen::VectorXd >& residuals )
{
    return std::make_shared< SingleObservationSet< double, double > >( observableType,
                                                                       linkDefinition,
                                                                       observations,
                                                                       times,
                                                                       receiver,
                                                                       std::vector< Eigen::VectorXd >( ),
                                                                       nullptr,
                                                                       nullptr,
                                                                       weights,
                                                                       residuals );
}

BOOST_AUTO_TEST_CASE( test_dataset_from_single_observation_set )
{
    const LinkDefinition linkDefinition = createOneWayLinkDefinition( "Station1" );

    std::vector< Eigen::VectorXd > observations;
    observations.push_back( ( Eigen::Vector2d( ) << 20.0, 21.0 ).finished( ) );
    observations.push_back( ( Eigen::Vector2d( ) << 22.0, 23.0 ).finished( ) );

    std::vector< Eigen::VectorXd > weights;
    weights.push_back( ( Eigen::Vector2d( ) << 4.0, 5.0 ).finished( ) );
    weights.push_back( ( Eigen::Vector2d( ) << 6.0, 7.0 ).finished( ) );

    std::vector< Eigen::VectorXd > residuals;
    residuals.push_back( ( Eigen::Vector2d( ) << 0.4, 0.5 ).finished( ) );
    residuals.push_back( ( Eigen::Vector2d( ) << 0.6, 0.7 ).finished( ) );

    std::shared_ptr< SingleObservationSet< double, double > > observationSet =
            createManualObservationSet( angular_position, linkDefinition, observations, { 3.0, 4.0 }, weights, residuals );
    std::shared_ptr< ObservationDataset< double, double > > dataset = createObservationDataset( observationSet );

    BOOST_CHECK_EQUAL( dataset->getNumberOfObservationSets( ), 1 );
    BOOST_CHECK_EQUAL( dataset->getNumberOfObservations( ), 2 );
    BOOST_CHECK_EQUAL( dataset->getTotalScalarSize( ), 4 );
    BOOST_CHECK_EQUAL( dataset->getNumberOfLinkDefinitions( ), 1 );
    BOOST_CHECK( dataset->getLinkDefinition( dataset->getObservationSetMetadata( 0 ).linkDefinitionId_ ) == linkDefinition );

    BOOST_CHECK_EQUAL( dataset->getObservationRow( 1 ).firstScalarComponent_, 2 );
    BOOST_CHECK_EQUAL( dataset->getObservationRow( 1 ).scalarSize_, 2 );
    BOOST_CHECK_EQUAL( dataset->getScalarComponentRow( 3 ).observationId_, 1 );
    BOOST_CHECK_EQUAL( dataset->getScalarComponentRow( 3 ).componentIndex_, 1 );

    const EstimationVectorProjection< double, double > projection = dataset->createLegacyProjection( );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( projection.getObservationVector( ), observationSet->getObservationsVector( ), 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( projection.getResidualVector( ), observationSet->getResidualsVector( ), 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( projection.getWeightVector( ), observationSet->getWeightsVector( ), 1.0E-15 );
    std::vector< double > legacyTimes;
    for( const double observationTime : observationSet->getObservationTimes( ) )
    {
        for( unsigned int i = 0; i < observationSet->getSingleObservableSize( ); ++i )
        {
            legacyTimes.push_back( observationTime );
        }
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(
            projection.getTimes( ).begin( ), projection.getTimes( ).end( ), legacyTimes.begin( ), legacyTimes.end( ) );
}

BOOST_AUTO_TEST_CASE( test_dataset_from_observation_collection_uses_legacy_order )
{
    const LinkDefinition station1LinkDefinition = createOneWayLinkDefinition( "Station1" );
    const LinkDefinition station2LinkDefinition = createOneWayLinkDefinition( "Station2" );

    std::shared_ptr< SingleObservationSet< double, double > > rangeSet =
            createManualObservationSet( one_way_range,
                                        station1LinkDefinition,
                                        { ( Eigen::VectorXd( 1 ) << 10.0 ).finished( ), ( Eigen::VectorXd( 1 ) << 11.0 ).finished( ) },
                                        { 1.0, 2.0 },
                                        { ( Eigen::VectorXd( 1 ) << 2.0 ).finished( ), ( Eigen::VectorXd( 1 ) << 3.0 ).finished( ) },
                                        { ( Eigen::VectorXd( 1 ) << 0.1 ).finished( ), ( Eigen::VectorXd( 1 ) << 0.2 ).finished( ) } );

    std::shared_ptr< SingleObservationSet< double, double > > angularSet = createManualObservationSet(
            angular_position,
            station1LinkDefinition,
            { ( Eigen::Vector2d( ) << 20.0, 21.0 ).finished( ), ( Eigen::Vector2d( ) << 22.0, 23.0 ).finished( ) },
            { 3.0, 4.0 },
            { ( Eigen::Vector2d( ) << 4.0, 5.0 ).finished( ), ( Eigen::Vector2d( ) << 6.0, 7.0 ).finished( ) },
            { ( Eigen::Vector2d( ) << 0.4, 0.5 ).finished( ), ( Eigen::Vector2d( ) << 0.6, 0.7 ).finished( ) } );

    std::shared_ptr< SingleObservationSet< double, double > > positionSet =
            createManualObservationSet( position_observable,
                                        station2LinkDefinition,
                                        { ( Eigen::Vector3d( ) << 30.0, 31.0, 32.0 ).finished( ) },
                                        { 5.0 },
                                        { ( Eigen::Vector3d( ) << 8.0, 9.0, 10.0 ).finished( ) },
                                        { ( Eigen::Vector3d( ) << 0.8, 0.9, 1.0 ).finished( ) } );

    std::shared_ptr< ObservationCollection< double, double > > observationCollection =
            std::make_shared< ObservationCollection< double, double > >(
                    std::vector< std::shared_ptr< SingleObservationSet< double, double > > >( { positionSet, angularSet, rangeSet } ) );

    std::shared_ptr< ObservationDataset< double, double > > dataset = createObservationDataset( observationCollection );
    const EstimationVectorProjection< double, double > projection = dataset->createLegacyProjection( );

    BOOST_CHECK_EQUAL( dataset->getNumberOfObservationSets( ), 3 );
    BOOST_CHECK_EQUAL( dataset->getNumberOfObservations( ), 5 );
    BOOST_CHECK_EQUAL( dataset->getTotalScalarSize( ), observationCollection->getTotalObservableSize( ) );
    BOOST_CHECK_EQUAL( dataset->getNumberOfLinkDefinitions( ), 2 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( projection.getObservationVector( ), observationCollection->getObservationVector( ), 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( projection.getResidualVector( ), observationCollection->getConcatenatedResiduals( ), 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( projection.getWeightVector( ), observationCollection->getConcatenatedWeights( ), 1.0E-15 );
    const std::vector< double > legacyTimes = observationCollection->getConcatenatedTimeVector( );
    BOOST_CHECK_EQUAL_COLLECTIONS(
            projection.getTimes( ).begin( ), projection.getTimes( ).end( ), legacyTimes.begin( ), legacyTimes.end( ) );

    BOOST_CHECK_EQUAL( dataset->getObservationSetMetadata( 0 ).observableType_, one_way_range );
    BOOST_CHECK_EQUAL( dataset->getObservationSetMetadata( 1 ).observableType_, angular_position );
    BOOST_CHECK_EQUAL( dataset->getObservationSetMetadata( 2 ).observableType_, position_observable );
    BOOST_CHECK_EQUAL( dataset->getObservationSetMetadata( 0 ).linkDefinitionId_,
                       dataset->getObservationSetMetadata( 1 ).linkDefinitionId_ );
}

BOOST_AUTO_TEST_CASE( test_legacy_collection_getters_use_dataset_backend_after_mutation )
{
    const LinkDefinition linkDefinition = createOneWayLinkDefinition( "Station1" );

    std::shared_ptr< SingleObservationSet< double, double > > rangeSet =
            createManualObservationSet( one_way_range,
                                        linkDefinition,
                                        { ( Eigen::VectorXd( 1 ) << 10.0 ).finished( ), ( Eigen::VectorXd( 1 ) << 11.0 ).finished( ) },
                                        { 1.0, 2.0 },
                                        { ( Eigen::VectorXd( 1 ) << 2.0 ).finished( ), ( Eigen::VectorXd( 1 ) << 3.0 ).finished( ) },
                                        { ( Eigen::VectorXd( 1 ) << 0.1 ).finished( ), ( Eigen::VectorXd( 1 ) << 0.2 ).finished( ) } );

    std::shared_ptr< SingleObservationSet< double, double > > angularSet = createManualObservationSet(
            angular_position,
            linkDefinition,
            { ( Eigen::Vector2d( ) << 20.0, 21.0 ).finished( ), ( Eigen::Vector2d( ) << 22.0, 23.0 ).finished( ) },
            { 3.0, 4.0 },
            { ( Eigen::Vector2d( ) << 4.0, 5.0 ).finished( ), ( Eigen::Vector2d( ) << 6.0, 7.0 ).finished( ) },
            { ( Eigen::Vector2d( ) << 0.4, 0.5 ).finished( ), ( Eigen::Vector2d( ) << 0.6, 0.7 ).finished( ) } );

    ObservationCollection< double, double > observationCollection(
            std::vector< std::shared_ptr< SingleObservationSet< double, double > > >( { angularSet, rangeSet } ) );

    Eigen::VectorXd expectedInitialObservationVector( 6 );
    expectedInitialObservationVector << 10.0, 11.0, 20.0, 21.0, 22.0, 23.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( observationCollection.getObservationVector( ), expectedInitialObservationVector, 1.0E-15 );

    Eigen::VectorXd replacementObservationVector( 6 );
    replacementObservationVector << 100.0, 101.0, 200.0, 201.0, 202.0, 203.0;
    observationCollection.setObservations( replacementObservationVector );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( observationCollection.getObservationVector( ), replacementObservationVector, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( observationCollection.getObservationVectorReference( ), replacementObservationVector, 1.0E-15 );

    const std::vector< double > expectedScalarTimes = { 1.0, 2.0, 3.0, 3.0, 4.0, 4.0 };
    const std::vector< double > legacyScalarTimes = observationCollection.getConcatenatedTimeVector( );
    BOOST_CHECK_EQUAL_COLLECTIONS(
            legacyScalarTimes.begin( ), legacyScalarTimes.end( ), expectedScalarTimes.begin( ), expectedScalarTimes.end( ) );
}

BOOST_AUTO_TEST_CASE( test_dataset_parser_selection )
{
    ObservationDataset< double, double > dataset;

    const LinkDefinition station1LinkDefinition = createOneWayLinkDefinition( "Station1" );
    const LinkDefinition station2LinkDefinition = createOneWayLinkDefinition( "Station2" );

    LinkEnds marsLinkEnds;
    marsLinkEnds[ transmitter ] = LinkEndId( "Mars", "Station3" );
    marsLinkEnds[ receiver ] = LinkEndId( "Vehicle", "" );
    const LinkDefinition marsLinkDefinition( marsLinkEnds );

    const ObservationSetId rangeSetId =
            dataset.addObservationSet( one_way_range,
                                       station1LinkDefinition,
                                       { ( Eigen::VectorXd( 1 ) << 10.0 ).finished( ), ( Eigen::VectorXd( 1 ) << 11.0 ).finished( ) },
                                       { 1.0, 2.0 },
                                       receiver );
    const ObservationSetId angularSetId = dataset.addObservationSet(
            angular_position,
            station2LinkDefinition,
            { ( Eigen::Vector2d( ) << 20.0, 21.0 ).finished( ), ( Eigen::Vector2d( ) << 22.0, 23.0 ).finished( ) },
            { 5.0, 6.0 },
            receiver );
    const ObservationSetId positionSetId = dataset.addObservationSet(
            position_observable, marsLinkDefinition, { ( Eigen::Vector3d( ) << 30.0, 31.0, 32.0 ).finished( ) }, { 9.0 }, receiver );

    std::vector< ObservationSetId > selectedSetIds = dataset.getObservationSetIds( observationParser( one_way_range ) );
    BOOST_CHECK_EQUAL_COLLECTIONS( selectedSetIds.begin( ), selectedSetIds.end( ), &rangeSetId, &rangeSetId + 1 );

    selectedSetIds = dataset.getObservationSetIds( observationParser( angular_position, true ) );
    const std::vector< ObservationSetId > nonAngularSetIds = { rangeSetId, positionSetId };
    BOOST_CHECK_EQUAL_COLLECTIONS( selectedSetIds.begin( ), selectedSetIds.end( ), nonAngularSetIds.begin( ), nonAngularSetIds.end( ) );

    selectedSetIds = dataset.getObservationSetIds( observationParser( station2LinkDefinition.linkEnds_ ) );
    BOOST_CHECK_EQUAL_COLLECTIONS( selectedSetIds.begin( ), selectedSetIds.end( ), &angularSetId, &angularSetId + 1 );

    selectedSetIds = dataset.getObservationSetIds( observationParser( std::string( "Station2" ), true ) );
    BOOST_CHECK_EQUAL_COLLECTIONS( selectedSetIds.begin( ), selectedSetIds.end( ), &angularSetId, &angularSetId + 1 );

    selectedSetIds = dataset.getObservationSetIds( observationParser( std::string( "Earth" ) ) );
    const std::vector< ObservationSetId > earthSetIds = { rangeSetId, angularSetId };
    BOOST_CHECK_EQUAL_COLLECTIONS( selectedSetIds.begin( ), selectedSetIds.end( ), earthSetIds.begin( ), earthSetIds.end( ) );

    selectedSetIds = dataset.getObservationSetIds( observationParser( std::make_pair( 0.0, 3.0 ) ) );
    BOOST_CHECK_EQUAL_COLLECTIONS( selectedSetIds.begin( ), selectedSetIds.end( ), &rangeSetId, &rangeSetId + 1 );

    selectedSetIds = dataset.getObservationSetIds(
            observationParser( std::vector< std::shared_ptr< ObservationCollectionParser > >(
                                       { observationParser( one_way_range ), observationParser( std::make_pair( 0.0, 3.0 ) ) } ),
                               true ) );
    BOOST_CHECK_EQUAL_COLLECTIONS( selectedSetIds.begin( ), selectedSetIds.end( ), &rangeSetId, &rangeSetId + 1 );

    selectedSetIds = dataset.getObservationSetIds(
            observationParser( std::vector< std::shared_ptr< ObservationCollectionParser > >(
                                       { observationParser( one_way_range ), observationParser( std::string( "Mars" ) ) } ),
                               false ) );
    const std::vector< ObservationSetId > unionSetIds = { rangeSetId, positionSetId };
    BOOST_CHECK_EQUAL_COLLECTIONS( selectedSetIds.begin( ), selectedSetIds.end( ), unionSetIds.begin( ), unionSetIds.end( ) );
}

BOOST_AUTO_TEST_CASE( test_dataset_empty_sets_and_invalid_inputs )
{
    const LinkDefinition linkDefinition = createOneWayLinkDefinition( "Station1" );
    ObservationDataset< double, double > dataset;

    const ObservationSetId emptySetId = dataset.addObservationSet(
            one_way_range, linkDefinition, std::vector< Eigen::VectorXd >( ), std::vector< double >( ), receiver );

    BOOST_CHECK_EQUAL( dataset.getNumberOfObservationSets( ), 1 );
    BOOST_CHECK_EQUAL( dataset.getNumberOfObservationsForSet( emptySetId ), 0 );
    BOOST_CHECK_EQUAL( dataset.getTotalScalarSizeForSet( emptySetId ), 0 );
    BOOST_CHECK_EQUAL( dataset.createLegacyProjection( ).getObservationVector( ).size( ), 0 );
    BOOST_CHECK( std::isnan( dataset.getTimeBoundsForSet( emptySetId ).first ) );
    BOOST_CHECK( std::isnan( dataset.getTimeBoundsForSet( emptySetId ).second ) );

    BOOST_CHECK_THROW(
            dataset.addObservationSet(
                    one_way_range, linkDefinition, { ( Eigen::VectorXd( 1 ) << 1.0 ).finished( ) }, std::vector< double >( ), receiver ),
            std::runtime_error );

    BOOST_CHECK_THROW( dataset.addObservationSet( one_way_range,
                                                  linkDefinition,
                                                  { ( Eigen::VectorXd( 1 ) << 1.0 ).finished( ) },
                                                  { 1.0 },
                                                  receiver,
                                                  std::vector< Eigen::VectorXd >( ),
                                                  nullptr,
                                                  nullptr,
                                                  { ( Eigen::Vector2d( ) << 1.0, 2.0 ).finished( ) } ),
                       std::runtime_error );

    BOOST_CHECK_THROW( dataset.setObservationVectorForSet( emptySetId, ( Eigen::VectorXd( 1 ) << 1.0 ).finished( ) ), std::runtime_error );
    BOOST_CHECK_THROW( dataset.removeObservationFromSet( emptySetId, 0 ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( test_dataset_duplicate_filter_and_move_edge_cases )
{
    const LinkDefinition linkDefinition = createOneWayLinkDefinition( "Station1" );

    ObservationDataset< double, double > dataset;
    const ObservationSetId sourceSetId = dataset.addObservationSet( one_way_range,
                                                                    linkDefinition,
                                                                    { ( Eigen::VectorXd( 1 ) << 10.0 ).finished( ),
                                                                      ( Eigen::VectorXd( 1 ) << 11.0 ).finished( ),
                                                                      ( Eigen::VectorXd( 1 ) << 11.0 ).finished( ),
                                                                      ( Eigen::VectorXd( 1 ) << 13.0 ).finished( ) },
                                                                    { 4.0, 1.0, 1.0, 3.0 },
                                                                    receiver,
                                                                    std::vector< Eigen::VectorXd >( ),
                                                                    nullptr,
                                                                    nullptr,
                                                                    { ( Eigen::VectorXd( 1 ) << 2.0 ).finished( ),
                                                                      ( Eigen::VectorXd( 1 ) << 3.0 ).finished( ),
                                                                      ( Eigen::VectorXd( 1 ) << 4.0 ).finished( ),
                                                                      ( Eigen::VectorXd( 1 ) << 5.0 ).finished( ) },
                                                                    { ( Eigen::VectorXd( 1 ) << 0.1 ).finished( ),
                                                                      ( Eigen::VectorXd( 1 ) << 0.2 ).finished( ),
                                                                      ( Eigen::VectorXd( 1 ) << 0.2 ).finished( ),
                                                                      ( Eigen::VectorXd( 1 ) << 0.3 ).finished( ) },
                                                                    true,
                                                                    true );

    BOOST_CHECK_EQUAL( dataset.getNumberOfObservationsForSet( sourceSetId ), 3 );
    const std::vector< double > sortedUniqueTimes = dataset.getObservationTimesForSet( sourceSetId );
    const std::vector< double > expectedSortedUniqueTimes = { 1.0, 3.0, 4.0 };
    BOOST_CHECK_EQUAL_COLLECTIONS(
            sortedUniqueTimes.begin( ), sortedUniqueTimes.end( ), expectedSortedUniqueTimes.begin( ), expectedSortedUniqueTimes.end( ) );

    const std::vector< unsigned int > residualFilteredIndices = dataset.getFilteredObservationIndices(
            sourceSetId, std::make_shared< ObservationFilter< double > >( residual_filtering, 0.25 ) );
    const std::vector< unsigned int > expectedResidualFilteredIndices = { 1 };
    BOOST_CHECK_EQUAL_COLLECTIONS( residualFilteredIndices.begin( ),
                                   residualFilteredIndices.end( ),
                                   expectedResidualFilteredIndices.begin( ),
                                   expectedResidualFilteredIndices.end( ) );

    const std::vector< unsigned int > absoluteValueFilteredIndices = dataset.getFilteredObservationIndices(
            sourceSetId, std::make_shared< ObservationFilter< double > >( absolute_value_filtering, 12.0 ) );
    const std::vector< unsigned int > expectedAbsoluteValueFilteredIndices = { 1 };
    BOOST_CHECK_EQUAL_COLLECTIONS( absoluteValueFilteredIndices.begin( ),
                                   absoluteValueFilteredIndices.end( ),
                                   expectedAbsoluteValueFilteredIndices.begin( ),
                                   expectedAbsoluteValueFilteredIndices.end( ) );

    ObservationDataset< double, double > targetDataset;
    const ObservationSetId targetSetId = targetDataset.addObservationSet(
            one_way_range, linkDefinition, std::vector< Eigen::VectorXd >( ), std::vector< double >( ), receiver );

    dataset.moveObservationsToSet( sourceSetId, targetDataset, targetSetId, std::vector< unsigned int >( { 0, 2 } ), false );
    BOOST_CHECK_EQUAL( dataset.getNumberOfObservationsForSet( sourceSetId ), 3 );
    BOOST_CHECK_EQUAL( targetDataset.getNumberOfObservationsForSet( targetSetId ), 2 );

    dataset.moveObservationsToSet( sourceSetId, targetDataset, targetSetId, std::vector< unsigned int >( { 1 } ), true );
    BOOST_CHECK_EQUAL( dataset.getNumberOfObservationsForSet( sourceSetId ), 2 );
    BOOST_CHECK_EQUAL( targetDataset.getNumberOfObservationsForSet( targetSetId ), 3 );

    const Eigen::VectorXd targetObservationVector = targetDataset.getObservationVectorForSet( targetSetId );
    Eigen::VectorXd expectedTargetObservationVector( 3 );
    expectedTargetObservationVector << 11.0, 13.0, 10.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( targetObservationVector, expectedTargetObservationVector, 1.0E-15 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
