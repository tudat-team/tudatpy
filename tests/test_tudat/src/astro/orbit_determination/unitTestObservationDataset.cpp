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
#include "tudat/astro/observation_models/observationSimulator.h"
#include "tudat/astro/orbit_determination/podInputOutputTypes.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManagerHelpers.h"
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

class ZeroObservationSimulator : public ObservationSimulatorBase< double, double >
{
public:
    ZeroObservationSimulator( const ObservableType observableType, const int observableSize ):
        ObservationSimulatorBase< double, double >( observableType ), observableSize_( observableSize )
    {}

    int getObservationSize( ) override
    {
        return observableSize_;
    }

    void computeObservations( const std::vector< double >& times,
                              const LinkEnds,
                              const LinkEndType,
                              const std::shared_ptr< ObservationAncillarySimulationSettings >,
                              Eigen::Matrix< double, Eigen::Dynamic, 1 >& observationsVector ) override
    {
        observationsVector = Eigen::VectorXd::Zero( static_cast< int >( times.size( ) ) * observableSize_ );
    }

private:
    int observableSize_;
};

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

BOOST_AUTO_TEST_CASE( test_direct_dataset_legacy_order_bookkeeping_and_residuals )
{
    const LinkDefinition stationALinkDefinition = createOneWayLinkDefinition( "StationA" );
    const LinkDefinition stationBLinkDefinition = createOneWayLinkDefinition( "StationB" );

    ObservationDataset< double, double > dataset;
    const ObservationSetId positionSetId = dataset.addObservationSet( position_observable,
                                                                      stationBLinkDefinition,
                                                                      { ( Eigen::Vector3d( ) << 30.0, 31.0, 32.0 ).finished( ) },
                                                                      { 5.0 },
                                                                      receiver,
                                                                      std::vector< Eigen::VectorXd >( ),
                                                                      nullptr,
                                                                      nullptr,
                                                                      { ( Eigen::Vector3d( ) << 8.0, 9.0, 10.0 ).finished( ) },
                                                                      { ( Eigen::Vector3d( ) << 0.8, 0.9, 1.0 ).finished( ) } );
    const ObservationSetId angularSetId = dataset.addObservationSet(
            angular_position,
            stationALinkDefinition,
            { ( Eigen::Vector2d( ) << 20.0, 21.0 ).finished( ), ( Eigen::Vector2d( ) << 22.0, 23.0 ).finished( ) },
            { 3.0, 4.0 },
            receiver,
            std::vector< Eigen::VectorXd >( ),
            nullptr,
            nullptr,
            { ( Eigen::Vector2d( ) << 4.0, 5.0 ).finished( ), ( Eigen::Vector2d( ) << 6.0, 7.0 ).finished( ) },
            { ( Eigen::Vector2d( ) << 0.4, 0.5 ).finished( ), ( Eigen::Vector2d( ) << 0.6, 0.7 ).finished( ) } );
    const ObservationSetId rangeSetId =
            dataset.addObservationSet( one_way_range,
                                       stationALinkDefinition,
                                       { ( Eigen::VectorXd( 1 ) << 10.0 ).finished( ), ( Eigen::VectorXd( 1 ) << 11.0 ).finished( ) },
                                       { 1.0, 2.0 },
                                       receiver,
                                       std::vector< Eigen::VectorXd >( ),
                                       nullptr,
                                       nullptr,
                                       { ( Eigen::VectorXd( 1 ) << 2.0 ).finished( ), ( Eigen::VectorXd( 1 ) << 3.0 ).finished( ) },
                                       { ( Eigen::VectorXd( 1 ) << 0.1 ).finished( ), ( Eigen::VectorXd( 1 ) << 0.2 ).finished( ) } );

    BOOST_CHECK_EQUAL( positionSetId, 0 );
    BOOST_CHECK_EQUAL( angularSetId, 1 );
    BOOST_CHECK_EQUAL( rangeSetId, 2 );

    const std::vector< ObservationSetId > expectedLegacySetIds = { rangeSetId, angularSetId, positionSetId };
    const std::vector< ObservationSetId > legacySetIds = dataset.getSetIdsInLegacyOrder( );
    BOOST_CHECK_EQUAL_COLLECTIONS( legacySetIds.begin( ), legacySetIds.end( ), expectedLegacySetIds.begin( ), expectedLegacySetIds.end( ) );

    Eigen::VectorXd expectedLegacyObservations( 9 );
    expectedLegacyObservations << 10.0, 11.0, 20.0, 21.0, 22.0, 23.0, 30.0, 31.0, 32.0;
    const EstimationVectorProjection< double, double > legacyProjection = dataset.createLegacyProjection( );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( legacyProjection.getObservationVector( ), expectedLegacyObservations, 1.0E-15 );

    const std::vector< double > expectedLegacyTimes = { 1.0, 2.0, 3.0, 3.0, 4.0, 4.0, 5.0, 5.0, 5.0 };
    const std::vector< double >& legacyTimes = legacyProjection.getTimes( );
    BOOST_CHECK_EQUAL_COLLECTIONS( legacyTimes.begin( ), legacyTimes.end( ), expectedLegacyTimes.begin( ), expectedLegacyTimes.end( ) );

    const std::vector< std::pair< int, int > > expectedLegacyStartAndSize = { { 0, 2 }, { 2, 4 }, { 6, 3 } };
    const std::vector< std::pair< int, int > > legacyStartAndSize = dataset.getObservationSetStartAndSize( );
    BOOST_REQUIRE_EQUAL( legacyStartAndSize.size( ), expectedLegacyStartAndSize.size( ) );
    for( unsigned int i = 0; i < legacyStartAndSize.size( ); ++i )
    {
        BOOST_CHECK_EQUAL( legacyStartAndSize.at( i ).first, expectedLegacyStartAndSize.at( i ).first );
        BOOST_CHECK_EQUAL( legacyStartAndSize.at( i ).second, expectedLegacyStartAndSize.at( i ).second );
    }

    const std::vector< std::pair< int, int > > expectedDatasetStartAndSize = { { 0, 3 }, { 3, 4 }, { 7, 2 } };
    const std::vector< std::pair< int, int > > datasetStartAndSize = dataset.getObservationSetStartAndSizeInDatasetOrder( );
    BOOST_REQUIRE_EQUAL( datasetStartAndSize.size( ), expectedDatasetStartAndSize.size( ) );
    for( unsigned int i = 0; i < datasetStartAndSize.size( ); ++i )
    {
        BOOST_CHECK_EQUAL( datasetStartAndSize.at( i ).first, expectedDatasetStartAndSize.at( i ).first );
        BOOST_CHECK_EQUAL( datasetStartAndSize.at( i ).second, expectedDatasetStartAndSize.at( i ).second );
    }

    Eigen::VectorXd expectedLegacyWeights( 9 );
    expectedLegacyWeights << 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( legacyProjection.getWeightVector( ), expectedLegacyWeights, 1.0E-15 );

    std::map< ObservableType, std::shared_ptr< ObservationSimulatorBase< double, double > > > simulators;
    simulators[ one_way_range ] = std::make_shared< ZeroObservationSimulator >( one_way_range, 1 );
    simulators[ angular_position ] = std::make_shared< ZeroObservationSimulator >( angular_position, 2 );
    simulators[ position_observable ] = std::make_shared< ZeroObservationSimulator >( position_observable, 3 );

    Eigen::VectorXd residuals;
    simulation_setup::calculateResiduals< double, double >(
            std::make_shared< ObservationDataset< double, double > >( dataset ), simulators, residuals );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( residuals, expectedLegacyObservations, 1.0E-15 );
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

BOOST_AUTO_TEST_CASE( test_covariance_input_tracks_legacy_collection_dataset_mutations )
{
    const LinkDefinition linkDefinition = createOneWayLinkDefinition( "Station1" );

    std::shared_ptr< SingleObservationSet< double, double > > rangeSet =
            createManualObservationSet( one_way_range,
                                        linkDefinition,
                                        { ( Eigen::VectorXd( 1 ) << 10.0 ).finished( ),
                                          ( Eigen::VectorXd( 1 ) << 11.0 ).finished( ),
                                          ( Eigen::VectorXd( 1 ) << 12.0 ).finished( ) },
                                        { 1.0, 2.0, 3.0 },
                                        { ( Eigen::VectorXd( 1 ) << 1.0 ).finished( ),
                                          ( Eigen::VectorXd( 1 ) << 1.0 ).finished( ),
                                          ( Eigen::VectorXd( 1 ) << 1.0 ).finished( ) },
                                        { ( Eigen::VectorXd( 1 ) << 0.0 ).finished( ),
                                          ( Eigen::VectorXd( 1 ) << 0.0 ).finished( ),
                                          ( Eigen::VectorXd( 1 ) << 0.0 ).finished( ) } );

    std::shared_ptr< ObservationCollection< double, double > > observationCollection =
            std::make_shared< ObservationCollection< double, double > >(
                    std::vector< std::shared_ptr< SingleObservationSet< double, double > > >( { rangeSet } ) );
    simulation_setup::CovarianceAnalysisInput< double, double > covarianceInput( observationCollection );

    std::shared_ptr< ObservationDataset< double, double > > inputDataset = covarianceInput.getObservationDataset( );
    BOOST_CHECK( inputDataset == observationCollection->getObservationDataset( ) );

    Eigen::VectorXd residuals( 3 );
    residuals << 0.1, 0.2, 0.3;
    observationCollection->setResiduals( residuals );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( inputDataset->createLegacyProjection( ).getResidualVector( ), residuals, 1.0E-15 );

    observationCollection->filterObservations(
            observationFilter( epochs_filtering, std::vector< double >( { 2.0 } ) ), observationParser( one_way_range ), false );

    BOOST_CHECK( inputDataset == observationCollection->getObservationDataset( ) );
    BOOST_CHECK_EQUAL( inputDataset->getTotalScalarSize( ), 2 );
    BOOST_CHECK_EQUAL( observationCollection->getTotalObservableSize( ), 2 );

    const std::vector< double > expectedTimes = { 1.0, 3.0 };
    const std::vector< double > inputTimes = inputDataset->createLegacyProjection( ).getTimes( );
    const std::vector< double > collectionTimes = observationCollection->getConcatenatedTimeVector( );
    BOOST_CHECK_EQUAL_COLLECTIONS( inputTimes.begin( ), inputTimes.end( ), expectedTimes.begin( ), expectedTimes.end( ) );
    BOOST_CHECK_EQUAL_COLLECTIONS( collectionTimes.begin( ), collectionTimes.end( ), expectedTimes.begin( ), expectedTimes.end( ) );
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

BOOST_AUTO_TEST_CASE( test_dataset_condition_viewer_rejection_and_reduced_dataset )
{
    const LinkDefinition linkDefinition = createOneWayLinkDefinition( "Station1" );
    ObservationDataset< double, double > dataset;

    dataset.addObservationSet( one_way_range,
                               linkDefinition,
                               { ( Eigen::VectorXd( 1 ) << 10.0 ).finished( ),
                                 ( Eigen::VectorXd( 1 ) << 20.0 ).finished( ),
                                 ( Eigen::VectorXd( 1 ) << 30.0 ).finished( ),
                                 ( Eigen::VectorXd( 1 ) << 40.0 ).finished( ) },
                               { 1.0, 2.0, 3.0, 4.0 },
                               receiver,
                               std::vector< Eigen::VectorXd >( ),
                               nullptr,
                               nullptr,
                               { ( Eigen::VectorXd( 1 ) << 1.0 ).finished( ),
                                 ( Eigen::VectorXd( 1 ) << 2.0 ).finished( ),
                                 ( Eigen::VectorXd( 1 ) << 3.0 ).finished( ),
                                 ( Eigen::VectorXd( 1 ) << 4.0 ).finished( ) },
                               { ( Eigen::VectorXd( 1 ) << 0.1 ).finished( ),
                                 ( Eigen::VectorXd( 1 ) << 0.2 ).finished( ),
                                 ( Eigen::VectorXd( 1 ) << 0.3 ).finished( ),
                                 ( Eigen::VectorXd( 1 ) << 0.4 ).finished( ) } );

    const ObservationCondition< double, double > middleTimes = ObservationCondition< double, double >::observableType( one_way_range ) &&
            ObservationCondition< double, double >::timeBounds( 1.5, 3.5 );
    const std::vector< ObservationId > selectedObservationIds = dataset.getObservationIdsMatchingCondition( middleTimes );
    const std::vector< ObservationId > expectedSelectedObservationIds = { 1, 2 };
    BOOST_CHECK_EQUAL_COLLECTIONS( selectedObservationIds.begin( ),
                                   selectedObservationIds.end( ),
                                   expectedSelectedObservationIds.begin( ),
                                   expectedSelectedObservationIds.end( ) );

    const ObservationDatasetViewer< double, double > viewer = dataset.createViewer( middleTimes );
    BOOST_CHECK_EQUAL( viewer.getNumberOfObservations( ), 2 );
    Eigen::VectorXd expectedViewerObservations( 2 );
    expectedViewerObservations << 20.0, 30.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( viewer.createEstimationProjection( ).getObservationVector( ), expectedViewerObservations, 1.0E-15 );

    const std::shared_ptr< ObservationDataset< double, double > > keptDataset = dataset.createNewAndKeep( middleTimes );
    BOOST_CHECK_EQUAL( keptDataset->getNumberOfObservationSets( ), 1 );
    BOOST_CHECK_EQUAL( keptDataset->getNumberOfObservations( ), 2 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            keptDataset->createLegacyProjection( ).getObservationVector( ), expectedViewerObservations, 1.0E-15 );

    const std::shared_ptr< ObservationDataset< double, double > > droppedDataset = dataset.createNewAndDrop( middleTimes );
    Eigen::VectorXd expectedDroppedObservations( 2 );
    expectedDroppedObservations << 10.0, 40.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            droppedDataset->createLegacyProjection( ).getObservationVector( ), expectedDroppedObservations, 1.0E-15 );

    dataset.rejectObservations( ObservationCondition< double, double >::timeBounds( 2.5, 3.5 ), "test rejection" );
    BOOST_CHECK( !dataset.getObservationRow( 2 ).isActive_ );
    BOOST_CHECK_EQUAL( dataset.getObservationRow( 2 ).rejectionReason_, "test rejection" );

    Eigen::VectorXd expectedActiveObservations( 3 );
    expectedActiveObservations << 10.0, 20.0, 40.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( dataset.createEstimationProjection( ).getObservationVector( ), expectedActiveObservations, 1.0E-15 );
    BOOST_CHECK_EQUAL( dataset.createLegacyProjection( ).getObservationVector( ).size( ), 4 );
    BOOST_CHECK_EQUAL( dataset.createViewer( ObservationCondition< double, double >::rejected( ) ).getNumberOfObservations( ), 1 );

    dataset.restoreObservations( ObservationCondition< double, double >::rejected( ) );
    BOOST_CHECK( dataset.getObservationRow( 2 ).isActive_ );
    BOOST_CHECK( dataset.getObservationRow( 2 ).rejectionReason_.empty( ) );
    BOOST_CHECK_EQUAL( dataset.createEstimationProjection( ).getObservationVector( ).size( ), 4 );

    const ObservationDatasetViewer< double, double > invalidatedViewer =
            dataset.createViewer( ObservationCondition< double, double >::all( ) );
    dataset.addObservationSet( one_way_range, linkDefinition, std::vector< Eigen::VectorXd >( ), std::vector< double >( ), receiver );
    BOOST_CHECK_THROW( invalidatedViewer.getNumberOfObservations( ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( test_dataset_compact_and_matrix_weights )
{
    const LinkDefinition linkDefinition = createOneWayLinkDefinition( "Station1" );

    ObservationDataset< double, double > scalarWeightDataset;
    const ObservationSetId scalarSetId = scalarWeightDataset.addObservationSetWithScalarWeight(
            angular_position,
            linkDefinition,
            { ( Eigen::Vector2d( ) << 10.0, 11.0 ).finished( ), ( Eigen::Vector2d( ) << 20.0, 21.0 ).finished( ) },
            { 1.0, 2.0 },
            receiver,
            5.0 );
    BOOST_CHECK_EQUAL( scalarSetId, 0 );
    Eigen::VectorXd expectedScalarWeights( 4 );
    expectedScalarWeights << 5.0, 5.0, 5.0, 5.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( scalarWeightDataset.getWeightVectorForSet( scalarSetId ), expectedScalarWeights, 1.0E-15 );
    const Eigen::MatrixXd expectedScalarWeightMatrix = 5.0 * Eigen::MatrixXd::Identity( 4, 4 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            scalarWeightDataset.createEstimationProjection( ).getWeightMatrix( ), expectedScalarWeightMatrix, 1.0E-15 );
    BOOST_CHECK( scalarWeightDataset.createEstimationProjection( ).isDiagonalWeightOnly( ) );

    Eigen::Matrix2d observationWeightBlock;
    observationWeightBlock << 2.0, 0.5, 0.5, 3.0;
    ObservationDataset< double, double > blockWeightDataset;
    const ObservationSetId blockSetId = blockWeightDataset.addObservationSetWithWeightBlock(
            angular_position,
            linkDefinition,
            { ( Eigen::Vector2d( ) << 10.0, 11.0 ).finished( ), ( Eigen::Vector2d( ) << 20.0, 21.0 ).finished( ) },
            { 1.0, 2.0 },
            receiver,
            observationWeightBlock );
    Eigen::MatrixXd expectedObservationBlockMatrix = Eigen::MatrixXd::Zero( 4, 4 );
    expectedObservationBlockMatrix.block( 0, 0, 2, 2 ) = observationWeightBlock;
    expectedObservationBlockMatrix.block( 2, 2, 2, 2 ) = observationWeightBlock;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            blockWeightDataset.createEstimationProjection( ).getWeightMatrix( ), expectedObservationBlockMatrix, 1.0E-15 );
    BOOST_CHECK( blockWeightDataset.createEstimationProjection( ).hasOffDiagonalWeights( ) );
    Eigen::VectorXd expectedBlockWeightDiagonal( 4 );
    expectedBlockWeightDiagonal << 2.0, 3.0, 2.0, 3.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( blockWeightDataset.getWeightVectorForSet( blockSetId ), expectedBlockWeightDiagonal, 1.0E-15 );

    Eigen::MatrixXd setWeightBlock = Eigen::MatrixXd::Zero( 4, 4 );
    setWeightBlock << 1.0, 0.1, 0.2, 0.3, 0.1, 2.0, 0.4, 0.5, 0.2, 0.4, 3.0, 0.6, 0.3, 0.5, 0.6, 4.0;
    ObservationDataset< double, double > setBlockWeightDataset;
    const ObservationSetId setBlockSetId = setBlockWeightDataset.addObservationSetWithSetWeightBlock(
            angular_position,
            linkDefinition,
            { ( Eigen::Vector2d( ) << 10.0, 11.0 ).finished( ), ( Eigen::Vector2d( ) << 20.0, 21.0 ).finished( ) },
            { 1.0, 2.0 },
            receiver,
            setWeightBlock );
    BOOST_CHECK( setBlockWeightDataset.hasWeightMatrixForSet( setBlockSetId ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( setBlockWeightDataset.getWeightMatrixForSet( setBlockSetId ), setWeightBlock, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( setBlockWeightDataset.createEstimationProjection( ).getWeightMatrix( ), setWeightBlock, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            setBlockWeightDataset.createEstimationProjection( ).getWeightVector( ), setWeightBlock.diagonal( ), 1.0E-15 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
