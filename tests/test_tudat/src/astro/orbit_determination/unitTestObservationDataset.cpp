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
#include "tudat/simulation/estimation_setup/observationOutputSettings.h"
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

void checkIds( const std::vector< ObservationSetId >& ids, const std::vector< ObservationSetId >& expectedIds )
{
    BOOST_CHECK_EQUAL_COLLECTIONS( ids.begin( ), ids.end( ), expectedIds.begin( ), expectedIds.end( ) );
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

BOOST_AUTO_TEST_CASE( test_dataset_storage_projection_and_residuals )
{
    const LinkDefinition stationALinkDefinition = createOneWayLinkDefinition( "StationA" );
    const LinkDefinition stationBLinkDefinition = createOneWayLinkDefinition( "StationB" );

    ObservationDataset< double, double > dataset;
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

    BOOST_CHECK_EQUAL( rangeSetId, 0 );
    BOOST_CHECK_EQUAL( angularSetId, 1 );
    BOOST_CHECK_EQUAL( positionSetId, 2 );
    BOOST_CHECK_EQUAL( dataset.getNumberOfObservationSets( ), 3 );
    BOOST_CHECK_EQUAL( dataset.getNumberOfObservations( ), 5 );
    BOOST_CHECK_EQUAL( dataset.getTotalScalarSize( ), 9 );
    BOOST_CHECK_EQUAL( dataset.getNumberOfLinkDefinitions( ), 2 );
    BOOST_CHECK_EQUAL( dataset.getObservationRow( 3 ).firstScalarComponent_, 4 );
    BOOST_CHECK_EQUAL( dataset.getObservationRow( 3 ).scalarSize_, 2 );
    BOOST_CHECK_EQUAL( dataset.getScalarComponentRow( 5 ).observationId_, 3 );
    BOOST_CHECK_EQUAL( dataset.getScalarComponentRow( 5 ).componentIndex_, 1 );

    Eigen::VectorXd expectedObservations( 9 );
    expectedObservations << 10.0, 11.0, 20.0, 21.0, 22.0, 23.0, 30.0, 31.0, 32.0;
    const EstimationVectorProjection< double, double > projection = dataset.createEstimationProjection( );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( projection.getObservationVector( ), expectedObservations, 1.0E-15 );

    const std::vector< double > expectedTimes = { 1.0, 2.0, 3.0, 3.0, 4.0, 4.0, 5.0, 5.0, 5.0 };
    const std::vector< double >& projectionTimes = projection.getTimes( );
    BOOST_CHECK_EQUAL_COLLECTIONS( projectionTimes.begin( ), projectionTimes.end( ), expectedTimes.begin( ), expectedTimes.end( ) );

    const std::vector< std::pair< int, int > > expectedDatasetStartAndSize = { { 0, 2 }, { 2, 4 }, { 6, 3 } };
    const std::vector< std::pair< int, int > > datasetStartAndSize = dataset.getObservationSetStartAndSizeInDatasetOrder( );
    BOOST_REQUIRE_EQUAL( datasetStartAndSize.size( ), expectedDatasetStartAndSize.size( ) );
    for( unsigned int i = 0; i < datasetStartAndSize.size( ); ++i )
    {
        BOOST_CHECK_EQUAL( datasetStartAndSize.at( i ).first, expectedDatasetStartAndSize.at( i ).first );
        BOOST_CHECK_EQUAL( datasetStartAndSize.at( i ).second, expectedDatasetStartAndSize.at( i ).second );
    }

    Eigen::VectorXd expectedWeights( 9 );
    expectedWeights << 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( projection.getWeightVector( ), expectedWeights, 1.0E-15 );

    std::map< ObservableType, std::shared_ptr< ObservationSimulatorBase< double, double > > > simulators;
    simulators[ one_way_range ] = std::make_shared< ZeroObservationSimulator >( one_way_range, 1 );
    simulators[ angular_position ] = std::make_shared< ZeroObservationSimulator >( angular_position, 2 );
    simulators[ position_observable ] = std::make_shared< ZeroObservationSimulator >( position_observable, 3 );

    Eigen::VectorXd residuals;
    simulation_setup::calculateResiduals< double, double >(
            std::make_shared< ObservationDataset< double, double > >( dataset ), simulators, residuals );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( residuals, expectedObservations, 1.0E-15 );
}

BOOST_AUTO_TEST_CASE( test_legacy_observation_interfaces_delegate_to_dataset_backend )
{
    const LinkDefinition station1LinkDefinition = createOneWayLinkDefinition( "Station1" );
    const LinkDefinition station2LinkDefinition = createOneWayLinkDefinition( "Station2" );

    const auto createManualObservationSet =
            []( const ObservableType observableType,
                const LinkDefinition& linkDefinition,
                const std::vector< Eigen::VectorXd >& observations,
                const std::vector< double >& times,
                const std::vector< Eigen::VectorXd >& weights,
                const std::vector< Eigen::VectorXd >& residuals ) -> std::shared_ptr< SingleObservationSet< double, double > > {
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
    };

    std::shared_ptr< SingleObservationSet< double, double > > angularConversionSet = createManualObservationSet(
            angular_position,
            station1LinkDefinition,
            { ( Eigen::Vector2d( ) << 20.0, 21.0 ).finished( ), ( Eigen::Vector2d( ) << 22.0, 23.0 ).finished( ) },
            { 3.0, 4.0 },
            { ( Eigen::Vector2d( ) << 4.0, 5.0 ).finished( ), ( Eigen::Vector2d( ) << 6.0, 7.0 ).finished( ) },
            { ( Eigen::Vector2d( ) << 0.4, 0.5 ).finished( ), ( Eigen::Vector2d( ) << 0.6, 0.7 ).finished( ) } );
    std::shared_ptr< ObservationDataset< double, double > > singleDataset = createObservationDataset( angularConversionSet );
    const EstimationVectorProjection< double, double > singleProjection = singleDataset->createLegacyProjection( );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( singleProjection.getObservationVector( ), angularConversionSet->getObservationsVector( ), 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( singleProjection.getResidualVector( ), angularConversionSet->getResidualsVector( ), 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( singleProjection.getWeightVector( ), angularConversionSet->getWeightsVector( ), 1.0E-15 );

    std::shared_ptr< SingleObservationSet< double, double > > singleSet =
            createManualObservationSet( one_way_range,
                                        station1LinkDefinition,
                                        { ( Eigen::VectorXd( 1 ) << 10.0 ).finished( ),
                                          ( Eigen::VectorXd( 1 ) << 11.0 ).finished( ),
                                          ( Eigen::VectorXd( 1 ) << 12.0 ).finished( ) },
                                        { 1.0, 2.0, 3.0 },
                                        { ( Eigen::VectorXd( 1 ) << 1.0 ).finished( ),
                                          ( Eigen::VectorXd( 1 ) << 1.0 ).finished( ),
                                          ( Eigen::VectorXd( 1 ) << 1.0 ).finished( ) },
                                        { ( Eigen::VectorXd( 1 ) << 0.1 ).finished( ),
                                          ( Eigen::VectorXd( 1 ) << 0.2 ).finished( ),
                                          ( Eigen::VectorXd( 1 ) << 0.3 ).finished( ) } );

    singleSet->filterObservations( observationFilter( epochs_filtering, std::vector< double >( { 2.0 } ) ), true );
    Eigen::VectorXd expectedSingleFilteredVector( 2 );
    expectedSingleFilteredVector << 10.0, 12.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( singleSet->getObservationsVector( ), expectedSingleFilteredVector, 1.0E-15 );
    Eigen::VectorXd expectedSavedFilteredVector( 1 );
    expectedSavedFilteredVector << 11.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            singleSet->getFilteredObservationSet( )->getObservationsVector( ), expectedSavedFilteredVector, 1.0E-15 );

    singleSet->filterObservations( observationFilter( epochs_filtering, std::vector< double >( { 2.0 } ), false ), true );
    Eigen::VectorXd expectedSingleRestoredVector( 3 );
    expectedSingleRestoredVector << 10.0, 11.0, 12.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( singleSet->getObservationsVector( ), expectedSingleRestoredVector, 1.0E-15 );

    singleSet->filterObservations( observationFilter( absolute_value_filtering, 11.0 ), false );
    Eigen::VectorXd expectedSingleDroppedVector( 2 );
    expectedSingleDroppedVector << 10.0, 11.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( singleSet->getObservationsVector( ), expectedSingleDroppedVector, 1.0E-15 );

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

    ObservationCollection< double, double > observationCollection(
            std::vector< std::shared_ptr< SingleObservationSet< double, double > > >( { positionSet, angularSet, rangeSet } ) );
    std::shared_ptr< ObservationDataset< double, double > > collectionDataset =
            createObservationDataset( std::make_shared< ObservationCollection< double, double > >( observationCollection ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            collectionDataset->createLegacyProjection( ).getObservationVector( ), observationCollection.getObservationVector( ), 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            collectionDataset->createLegacyProjection( ).getResidualVector( ), observationCollection.getConcatenatedResiduals( ), 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            collectionDataset->createLegacyProjection( ).getWeightVector( ), observationCollection.getConcatenatedWeights( ), 1.0E-15 );

    Eigen::VectorXd expectedInitialObservationVector( 9 );
    expectedInitialObservationVector << 10.0, 11.0, 20.0, 21.0, 22.0, 23.0, 30.0, 31.0, 32.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( observationCollection.getObservationVector( ), expectedInitialObservationVector, 1.0E-15 );

    const std::vector< ObservableType > observableTypes = observationCollection.getObservableTypes( );
    BOOST_CHECK_EQUAL( observableTypes.size( ), 3 );
    BOOST_CHECK_EQUAL( std::count( observableTypes.begin( ), observableTypes.end( ), one_way_range ), 1 );
    BOOST_CHECK_EQUAL( std::count( observableTypes.begin( ), observableTypes.end( ), angular_position ), 1 );
    BOOST_CHECK_EQUAL( std::count( observableTypes.begin( ), observableTypes.end( ), position_observable ), 1 );

    const std::vector< LinkEnds > collectionLinkEnds = observationCollection.getLinkEnds( );
    BOOST_CHECK_EQUAL( collectionLinkEnds.size( ), 2 );
    BOOST_CHECK_EQUAL( std::count( collectionLinkEnds.begin( ), collectionLinkEnds.end( ), station1LinkDefinition.linkEnds_ ), 1 );
    BOOST_CHECK_EQUAL( std::count( collectionLinkEnds.begin( ), collectionLinkEnds.end( ), station2LinkDefinition.linkEnds_ ), 1 );

    const std::vector< std::string > bodiesInLinkEnds = observationCollection.getBodiesInLinkEnds( );
    BOOST_CHECK_EQUAL( bodiesInLinkEnds.size( ), 2 );
    BOOST_CHECK_EQUAL( std::count( bodiesInLinkEnds.begin( ), bodiesInLinkEnds.end( ), "Earth" ), 1 );
    BOOST_CHECK_EQUAL( std::count( bodiesInLinkEnds.begin( ), bodiesInLinkEnds.end( ), "Vehicle" ), 1 );

    const std::vector< std::string > referencePointsInLinkEnds = observationCollection.getReferencePointsInLinkEnds( );
    BOOST_CHECK_EQUAL( referencePointsInLinkEnds.size( ), 2 );
    BOOST_CHECK_EQUAL( std::count( referencePointsInLinkEnds.begin( ), referencePointsInLinkEnds.end( ), "Station1" ), 1 );
    BOOST_CHECK_EQUAL( std::count( referencePointsInLinkEnds.begin( ), referencePointsInLinkEnds.end( ), "Station2" ), 1 );

    const std::vector< std::pair< double, double > > timeBoundsList = observationCollection.getTimeBoundsList( );
    BOOST_CHECK_EQUAL( timeBoundsList.size( ), 3 );
    BOOST_CHECK_EQUAL( std::count( timeBoundsList.begin( ), timeBoundsList.end( ), std::make_pair( 1.0, 2.0 ) ), 1 );
    BOOST_CHECK_EQUAL( std::count( timeBoundsList.begin( ), timeBoundsList.end( ), std::make_pair( 3.0, 4.0 ) ), 1 );
    BOOST_CHECK_EQUAL( std::count( timeBoundsList.begin( ), timeBoundsList.end( ), std::make_pair( 5.0, 5.0 ) ), 1 );

    std::vector< std::shared_ptr< SingleObservationSet< double, double > > > rangeObservationSets =
            observationCollection.getSingleObservationSets( observationParser( one_way_range ) );
    BOOST_REQUIRE_EQUAL( rangeObservationSets.size( ), 1 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            rangeObservationSets.at( 0 )->getObservationsVector( ), rangeSet->getObservationsVector( ), 1.0E-15 );

    std::vector< std::shared_ptr< SingleObservationSet< double, double > > > station1ObservationSets =
            observationCollection.getSingleObservationSets( observationParser( "Station1", true ) );
    BOOST_CHECK_EQUAL( station1ObservationSets.size( ), 2 );

    std::vector< std::shared_ptr< SingleObservationSet< double, double > > > station2ObservationSets =
            observationCollection.getSingleObservationSets( observationParser( station2LinkDefinition.linkEnds_ ) );
    BOOST_REQUIRE_EQUAL( station2ObservationSets.size( ), 1 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            station2ObservationSets.at( 0 )->getObservationsVector( ), positionSet->getObservationsVector( ), 1.0E-15 );

    std::vector< std::shared_ptr< ObservationCollectionParser > > station1RangeParserList;
    station1RangeParserList.push_back( observationParser( one_way_range ) );
    station1RangeParserList.push_back( observationParser( "Station1", true ) );
    std::vector< std::shared_ptr< SingleObservationSet< double, double > > > station1RangeObservationSets =
            observationCollection.getSingleObservationSets( observationParser( station1RangeParserList, true ) );
    BOOST_REQUIRE_EQUAL( station1RangeObservationSets.size( ), 1 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            station1RangeObservationSets.at( 0 )->getObservationsVector( ), rangeSet->getObservationsVector( ), 1.0E-15 );

    const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > rangeObservations =
            observationCollection.getObservations( observationParser( one_way_range ) );
    BOOST_REQUIRE_EQUAL( rangeObservations.size( ), 1 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( rangeObservations.at( 0 ), rangeSet->getObservationsVector( ), 1.0E-15 );

    const std::vector< std::vector< double > > rangeObservationTimes =
            observationCollection.getObservationTimes( observationParser( one_way_range ) );
    BOOST_REQUIRE_EQUAL( rangeObservationTimes.size( ), 1 );
    const std::vector< double > expectedRangeObservationTimes = rangeSet->getObservationTimes( );
    BOOST_CHECK_EQUAL_COLLECTIONS( rangeObservationTimes.at( 0 ).begin( ),
                                   rangeObservationTimes.at( 0 ).end( ),
                                   expectedRangeObservationTimes.begin( ),
                                   expectedRangeObservationTimes.end( ) );

    const std::map< ObservableType, std::pair< int, int > > observationTypeStartSize =
            observationCollection.getObservationTypeStartAndSize( );
    BOOST_CHECK_EQUAL( observationTypeStartSize.at( one_way_range ).first, 0 );
    BOOST_CHECK_EQUAL( observationTypeStartSize.at( one_way_range ).second, 2 );
    BOOST_CHECK_EQUAL( observationTypeStartSize.at( angular_position ).first, 2 );
    BOOST_CHECK_EQUAL( observationTypeStartSize.at( angular_position ).second, 4 );
    BOOST_CHECK_EQUAL( observationTypeStartSize.at( position_observable ).first, 6 );
    BOOST_CHECK_EQUAL( observationTypeStartSize.at( position_observable ).second, 3 );

    const std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > > observationSetStartSize =
            observationCollection.getObservationSetStartAndSize( );
    BOOST_CHECK_EQUAL( observationSetStartSize.at( one_way_range ).at( station1LinkDefinition.linkEnds_ ).at( 0 ).first, 0 );
    BOOST_CHECK_EQUAL( observationSetStartSize.at( one_way_range ).at( station1LinkDefinition.linkEnds_ ).at( 0 ).second, 2 );
    BOOST_CHECK_EQUAL( observationSetStartSize.at( angular_position ).at( station1LinkDefinition.linkEnds_ ).at( 0 ).first, 2 );
    BOOST_CHECK_EQUAL( observationSetStartSize.at( angular_position ).at( station1LinkDefinition.linkEnds_ ).at( 0 ).second, 4 );
    BOOST_CHECK_EQUAL( observationSetStartSize.at( position_observable ).at( station2LinkDefinition.linkEnds_ ).at( 0 ).first, 6 );
    BOOST_CHECK_EQUAL( observationSetStartSize.at( position_observable ).at( station2LinkDefinition.linkEnds_ ).at( 0 ).second, 3 );

    const std::vector< std::pair< int, int > > concatenatedStartSize = observationCollection.getConcatenatedObservationSetStartAndSize( );
    BOOST_REQUIRE_EQUAL( concatenatedStartSize.size( ), 3 );
    BOOST_CHECK_EQUAL( concatenatedStartSize.at( 0 ).first, 0 );
    BOOST_CHECK_EQUAL( concatenatedStartSize.at( 0 ).second, 2 );
    BOOST_CHECK_EQUAL( concatenatedStartSize.at( 1 ).first, 2 );
    BOOST_CHECK_EQUAL( concatenatedStartSize.at( 1 ).second, 4 );
    BOOST_CHECK_EQUAL( concatenatedStartSize.at( 2 ).first, 6 );
    BOOST_CHECK_EQUAL( concatenatedStartSize.at( 2 ).second, 3 );

    const std::vector< double > expectedEventTimes = { 1.0, 2.0, 3.0, 4.0, 5.0 };
    const std::vector< double > legacyEventTimes = observationCollection.getConcatenatedObservationTimes( );
    BOOST_CHECK_EQUAL_COLLECTIONS(
            legacyEventTimes.begin( ), legacyEventTimes.end( ), expectedEventTimes.begin( ), expectedEventTimes.end( ) );

    observationCollection.setConstantWeight( 9.0, observationParser( one_way_range ) );
    Eigen::VectorXd expectedWeightsAfterRangeConstant( 9 );
    expectedWeightsAfterRangeConstant << 9.0, 9.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( observationCollection.getConcatenatedWeights( ), expectedWeightsAfterRangeConstant, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( observationCollection.getObservationDataset( )->createLegacyProjection( ).getWeightVector( ),
                                       expectedWeightsAfterRangeConstant,
                                       1.0E-15 );

    Eigen::VectorXd replacementObservationVector( 9 );
    replacementObservationVector << 100.0, 101.0, 200.0, 201.0, 202.0, 203.0, 300.0, 301.0, 302.0;
    observationCollection.setObservations( replacementObservationVector );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( observationCollection.getObservationVector( ), replacementObservationVector, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( observationCollection.getObservationVectorReference( ), replacementObservationVector, 1.0E-15 );

    const std::vector< double > expectedScalarTimes = { 1.0, 2.0, 3.0, 3.0, 4.0, 4.0, 5.0, 5.0, 5.0 };
    const std::vector< double > legacyScalarTimes = observationCollection.getConcatenatedTimeVector( );
    BOOST_CHECK_EQUAL_COLLECTIONS(
            legacyScalarTimes.begin( ), legacyScalarTimes.end( ), expectedScalarTimes.begin( ), expectedScalarTimes.end( ) );

    observationCollection.filterObservations(
            observationFilter( epochs_filtering, std::vector< double >( { 2.0 } ) ), observationParser( one_way_range ), true );
    Eigen::VectorXd expectedCollectionFilteredVector( 8 );
    expectedCollectionFilteredVector << 100.0, 200.0, 201.0, 202.0, 203.0, 300.0, 301.0, 302.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( observationCollection.getObservationVector( ), expectedCollectionFilteredVector, 1.0E-15 );

    observationCollection.filterObservations(
            observationFilter( epochs_filtering, std::vector< double >( { 2.0 } ), false ), observationParser( one_way_range ), true );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( observationCollection.getObservationVector( ), replacementObservationVector, 1.0E-15 );

    std::shared_ptr< SingleObservationSet< double, double > > estimationRangeSet =
            createManualObservationSet( one_way_range,
                                        station1LinkDefinition,
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

    std::shared_ptr< ObservationCollection< double, double > > estimationCollection =
            std::make_shared< ObservationCollection< double, double > >(
                    std::vector< std::shared_ptr< SingleObservationSet< double, double > > >( { estimationRangeSet } ) );
    simulation_setup::CovarianceAnalysisInput< double, double > covarianceInput( estimationCollection );
    simulation_setup::EstimationInput< double, double > estimationInput( estimationCollection );

    std::shared_ptr< ObservationDataset< double, double > > covarianceInputDataset = covarianceInput.getObservationDataset( );
    std::shared_ptr< ObservationDataset< double, double > > estimationInputDataset = estimationInput.getObservationDataset( );
    BOOST_CHECK( covarianceInputDataset == estimationCollection->getObservationDataset( ) );
    BOOST_CHECK( estimationInputDataset == estimationCollection->getObservationDataset( ) );

    Eigen::VectorXd residuals( 3 );
    residuals << 0.1, 0.2, 0.3;
    estimationCollection->setResiduals( residuals );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( covarianceInputDataset->createLegacyProjection( ).getResidualVector( ), residuals, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( estimationInputDataset->createLegacyProjection( ).getResidualVector( ), residuals, 1.0E-15 );

    estimationCollection->filterObservations(
            observationFilter( epochs_filtering, std::vector< double >( { 2.0 } ) ), observationParser( one_way_range ), true );

    BOOST_CHECK( covarianceInputDataset == estimationCollection->getObservationDataset( ) );
    BOOST_CHECK( estimationInputDataset == estimationCollection->getObservationDataset( ) );
    BOOST_CHECK_EQUAL( covarianceInputDataset->getTotalScalarSize( ), 2 );
    BOOST_CHECK_EQUAL( estimationInputDataset->getTotalScalarSize( ), 2 );
    BOOST_CHECK_EQUAL( estimationCollection->getTotalObservableSize( ), 2 );

    const std::vector< double > expectedTimes = { 1.0, 3.0 };
    const std::vector< double > covarianceInputTimes = covarianceInputDataset->createLegacyProjection( ).getTimes( );
    const std::vector< double > estimationInputTimes = estimationInputDataset->createLegacyProjection( ).getTimes( );
    const std::vector< double > collectionTimes = estimationCollection->getConcatenatedTimeVector( );
    BOOST_CHECK_EQUAL_COLLECTIONS(
            covarianceInputTimes.begin( ), covarianceInputTimes.end( ), expectedTimes.begin( ), expectedTimes.end( ) );
    BOOST_CHECK_EQUAL_COLLECTIONS(
            estimationInputTimes.begin( ), estimationInputTimes.end( ), expectedTimes.begin( ), expectedTimes.end( ) );
    BOOST_CHECK_EQUAL_COLLECTIONS( collectionTimes.begin( ), collectionTimes.end( ), expectedTimes.begin( ), expectedTimes.end( ) );

    estimationCollection->filterObservations(
            observationFilter( epochs_filtering, std::vector< double >( { 2.0 } ), false ), observationParser( one_way_range ), true );

    BOOST_CHECK( covarianceInputDataset == estimationCollection->getObservationDataset( ) );
    BOOST_CHECK( estimationInputDataset == estimationCollection->getObservationDataset( ) );
    BOOST_CHECK_EQUAL( covarianceInputDataset->getTotalScalarSize( ), 3 );
    BOOST_CHECK_EQUAL( estimationInputDataset->getTotalScalarSize( ), 3 );
    const std::vector< double > expectedRestoredTimes = { 1.0, 2.0, 3.0 };
    const std::vector< double > restoredCovarianceInputTimes = covarianceInputDataset->createLegacyProjection( ).getTimes( );
    const std::vector< double > restoredEstimationInputTimes = estimationInputDataset->createLegacyProjection( ).getTimes( );
    BOOST_CHECK_EQUAL_COLLECTIONS( restoredCovarianceInputTimes.begin( ),
                                   restoredCovarianceInputTimes.end( ),
                                   expectedRestoredTimes.begin( ),
                                   expectedRestoredTimes.end( ) );
    BOOST_CHECK_EQUAL_COLLECTIONS( restoredEstimationInputTimes.begin( ),
                                   restoredEstimationInputTimes.end( ),
                                   expectedRestoredTimes.begin( ),
                                   expectedRestoredTimes.end( ) );
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
    BOOST_CHECK_EQUAL( dataset.createEstimationProjection( ).getObservationVector( ).size( ), 0 );
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

BOOST_AUTO_TEST_CASE( test_dataset_duplicate_selection_and_move_edge_cases )
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

    checkIds( dataset.getObservationIdsMatchingCondition( ObservationCondition< double, double >::residualAbsoluteValueGreaterThan(
                      ( Eigen::VectorXd( 1 ) << 0.25 ).finished( ) ) ),
              { 1 } );
    checkIds( dataset.getObservationIdsMatchingCondition( ObservationCondition< double, double >::observationAbsoluteValueGreaterThan(
                      ( Eigen::VectorXd( 1 ) << 12.0 ).finished( ) ) ),
              { 1 } );

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

BOOST_AUTO_TEST_CASE( test_dataset_row_conditions_cover_links_values_status_and_dependent_variables )
{
    const LinkDefinition station1LinkDefinition = createOneWayLinkDefinition( "Station1" );
    const LinkDefinition station2LinkDefinition = createOneWayLinkDefinition( "Station2" );

    std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > elevationSettings =
            simulation_setup::elevationAngleDependentVariable( receiver, LinkEndId( "Vehicle", "" ) );
    std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping > dependentVariableBookkeeping =
            std::make_shared< simulation_setup::ObservationDependentVariableBookkeeping >( one_way_range, station1LinkDefinition );
    dependentVariableBookkeeping->addDependentVariable( elevationSettings );

    ObservationDataset< double, double > dataset;
    dataset.addObservationSet( one_way_range,
                               station1LinkDefinition,
                               { ( Eigen::VectorXd( 1 ) << 10.0 ).finished( ),
                                 ( Eigen::VectorXd( 1 ) << 20.0 ).finished( ),
                                 ( Eigen::VectorXd( 1 ) << 30.0 ).finished( ) },
                               { 1.0, 2.0, 3.0 },
                               receiver,
                               { ( Eigen::VectorXd( 1 ) << 0.1 ).finished( ),
                                 ( Eigen::VectorXd( 1 ) << 0.4 ).finished( ),
                                 ( Eigen::VectorXd( 1 ) << 0.2 ).finished( ) },
                               dependentVariableBookkeeping,
                               nullptr,
                               std::vector< Eigen::VectorXd >( ),
                               { ( Eigen::VectorXd( 1 ) << 0.1 ).finished( ),
                                 ( Eigen::VectorXd( 1 ) << 0.6 ).finished( ),
                                 ( Eigen::VectorXd( 1 ) << 0.2 ).finished( ) } );

    dataset.addObservationSet( angular_position,
                               station2LinkDefinition,
                               { ( Eigen::Vector2d( ) << 1.0, 2.0 ).finished( ), ( Eigen::Vector2d( ) << 3.0, 4.0 ).finished( ) },
                               { 4.0, 5.0 },
                               receiver,
                               std::vector< Eigen::VectorXd >( ),
                               nullptr,
                               nullptr,
                               std::vector< Eigen::VectorXd >( ),
                               { ( Eigen::Vector2d( ) << 0.1, 0.2 ).finished( ), ( Eigen::Vector2d( ) << 0.7, 0.1 ).finished( ) } );

    checkIds(
            dataset.getObservationIdsMatchingCondition( ObservationCondition< double, double >::linkDefinition( station2LinkDefinition ) ),
            { 3, 4 } );
    checkIds( dataset.getObservationIdsMatchingCondition( ObservationCondition< double, double >::linkEndType( transmitter ) ),
              { 0, 1, 2, 3, 4 } );
    checkIds( dataset.getObservationIdsMatchingCondition(
                      ObservationCondition< double, double >::linkEnd( transmitter, LinkEndId( "Earth", "Station1" ) ) ),
              { 0, 1, 2 } );

    const ObservationCondition< double, double > rangeOrLastAngular =
            ObservationCondition< double, double >::observableType( one_way_range ) ||
            ObservationCondition< double, double >::timeBounds( 4.5, 5.5 );
    checkIds( dataset.getObservationIdsMatchingCondition( rangeOrLastAngular ), { 0, 1, 2, 4 } );

    checkIds( dataset.getObservationIdsMatchingCondition( !ObservationCondition< double, double >::timeBounds( 2.0, 4.0 ) ), { 0, 4 } );

    checkIds( dataset.getObservationIdsMatchingCondition( ObservationCondition< double, double >::observableType( one_way_range ) &&
                                                          ObservationCondition< double, double >::residualAbsoluteValueGreaterThan(
                                                                  ( Eigen::VectorXd( 1 ) << 0.5 ).finished( ) ) ),
              { 1 } );

    checkIds( dataset.getObservationIdsMatchingCondition( ObservationCondition< double, double >::observableType( angular_position ) &&
                                                          ObservationCondition< double, double >::observationAbsoluteValueGreaterThan(
                                                                  ( Eigen::Vector2d( ) << 2.5, 3.5 ).finished( ) ) ),
              { 4 } );

    checkIds( dataset.getObservationIdsMatchingCondition( ObservationCondition< double, double >::dependentVariableGreaterThan(
                      elevationSettings, ( Eigen::VectorXd( 1 ) << 0.3 ).finished( ) ) ),
              { 1 } );

    BOOST_CHECK_THROW( dataset.getObservationIdsMatchingCondition( ObservationCondition< double, double >::residualAbsoluteValueGreaterThan(
                               ( Eigen::Vector2d( ) << 0.1, 0.2 ).finished( ) ) ),
                       std::runtime_error );
    BOOST_CHECK_THROW( dataset.getObservationIdsMatchingCondition( ObservationCondition< double, double >::dependentVariableGreaterThan(
                               elevationSettings, ( Eigen::Vector2d( ) << 0.1, 0.2 ).finished( ) ) ),
                       std::runtime_error );

    dataset.rejectObservations(
            ObservationCondition< double, double >::observableType( one_way_range ) &&
                    ObservationCondition< double, double >::residualAbsoluteValueGreaterThan( ( Eigen::VectorXd( 1 ) << 0.5 ).finished( ) ),
            "large residual" );
    checkIds( dataset.getObservationIdsMatchingCondition( ObservationCondition< double, double >::rejected( ) ), { 1 } );
    checkIds( dataset.getObservationIdsMatchingCondition( ObservationCondition< double, double >::active( ) ), { 0, 2, 3, 4 } );
    dataset.restoreObservations( ObservationCondition< double, double >::rejected( ) );
    checkIds( dataset.getObservationIdsMatchingCondition( ObservationCondition< double, double >::active( ) ), { 0, 1, 2, 3, 4 } );
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
            keptDataset->createEstimationProjection( ).getObservationVector( ), expectedViewerObservations, 1.0E-15 );

    const std::shared_ptr< ObservationDataset< double, double > > droppedDataset = dataset.createNewAndDrop( middleTimes );
    Eigen::VectorXd expectedDroppedObservations( 2 );
    expectedDroppedObservations << 10.0, 40.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            droppedDataset->createEstimationProjection( ).getObservationVector( ), expectedDroppedObservations, 1.0E-15 );

    dataset.rejectObservations( ObservationCondition< double, double >::timeBounds( 2.5, 3.5 ), "test rejection" );
    BOOST_CHECK( !dataset.getObservationRow( 2 ).isActive_ );
    BOOST_CHECK_EQUAL( dataset.getObservationRow( 2 ).rejectionReason_, "test rejection" );

    Eigen::VectorXd expectedActiveObservations( 3 );
    expectedActiveObservations << 10.0, 20.0, 40.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( dataset.createEstimationProjection( ).getObservationVector( ), expectedActiveObservations, 1.0E-15 );
    BOOST_CHECK_EQUAL( dataset.createEstimationProjection( true ).getObservationVector( ).size( ), 4 );
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
