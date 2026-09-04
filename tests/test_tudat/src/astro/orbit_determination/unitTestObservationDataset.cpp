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
#include <iostream>
#include <sstream>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/observation_models/observationSimulator.h"
#include "tudat/astro/orbit_determination/podInputOutputTypes.h"
#include "tudat/math/basic/leastSquaresEstimation.h"
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

void checkIds( const std::vector< unsigned int >& ids, const std::vector< unsigned int >& expectedIds )
{
    // Selected observation ids must match exactly, including their dataset order.
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

/*!
 * Verifies the primary ObservationDataset storage model and flattened data output.
 *
 * Test outline: creates range, angular-position and position sets with explicit
 * weights and residuals. It checks set ids, row/scalar-component bookkeeping,
 * estimator-vector flattened data ordering, time replication for vector
 * observables, weight vector assembly and residual calculation through
 * observation simulators.
 */
BOOST_AUTO_TEST_CASE( test_dataset_storage_flattened_data_and_residuals )
{
    const LinkDefinition stationALinkDefinition = createOneWayLinkDefinition( "StationA" );
    const LinkDefinition stationBLinkDefinition = createOneWayLinkDefinition( "StationB" );

    ObservationDataset< double, double > dataset;
    const int rangeSetId =
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
    const int angularSetId = dataset.addObservationSet(
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
    const int positionSetId = dataset.addObservationSet( position_observable,
                                                         stationBLinkDefinition,
                                                         { ( Eigen::Vector3d( ) << 30.0, 31.0, 32.0 ).finished( ) },
                                                         { 5.0 },
                                                         receiver,
                                                         std::vector< Eigen::VectorXd >( ),
                                                         nullptr,
                                                         nullptr,
                                                         { ( Eigen::Vector3d( ) << 8.0, 9.0, 10.0 ).finished( ) },
                                                         { ( Eigen::Vector3d( ) << 0.8, 0.9, 1.0 ).finished( ) } );

    // Dataset-level counters and row bookkeeping must account for scalar and vector observables consistently.
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
    const FlattenedObservationData< double, double > flattenedData = dataset.createEstimationFlattenedObservationData( );

    // The estimator flattened data must concatenate scalar components in observation-set insertion order.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( flattenedData.getObservationVector( ), expectedObservations, 1.0E-15 );

    const std::vector< double > expectedTimes = { 1.0, 2.0, 3.0, 3.0, 4.0, 4.0, 5.0, 5.0, 5.0 };
    const std::vector< double >& flattenedDataTimes = flattenedData.getTimes( );

    // Vector-valued observables must repeat their event time once per scalar component.
    BOOST_CHECK_EQUAL_COLLECTIONS( flattenedDataTimes.begin( ), flattenedDataTimes.end( ), expectedTimes.begin( ), expectedTimes.end( ) );

    const std::vector< std::pair< int, int > > expectedDatasetStartAndSize = { { 0, 2 }, { 2, 4 }, { 6, 3 } };
    const std::vector< std::pair< int, int > > datasetStartAndSize = dataset.getObservationSetStartAndSizeInDatasetOrder( );

    // Set start/size data must describe the same scalar layout used by the estimator flattened data.
    BOOST_REQUIRE_EQUAL( datasetStartAndSize.size( ), expectedDatasetStartAndSize.size( ) );
    for( unsigned int i = 0; i < datasetStartAndSize.size( ); ++i )
    {
        BOOST_CHECK_EQUAL( datasetStartAndSize.at( i ).first, expectedDatasetStartAndSize.at( i ).first );
        BOOST_CHECK_EQUAL( datasetStartAndSize.at( i ).second, expectedDatasetStartAndSize.at( i ).second );
    }

    Eigen::VectorXd expectedWeights( 9 );
    expectedWeights << 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0;

    // Compact weights must be assembled in the same scalar-component order as the observations.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( flattenedData.getWeightVector( ), expectedWeights, 1.0E-15 );

    std::map< ObservableType, std::shared_ptr< ObservationSimulatorBase< double, double > > > simulators;
    simulators[ one_way_range ] = std::make_shared< ZeroObservationSimulator >( one_way_range, 1 );
    simulators[ angular_position ] = std::make_shared< ZeroObservationSimulator >( angular_position, 2 );
    simulators[ position_observable ] = std::make_shared< ZeroObservationSimulator >( position_observable, 3 );

    Eigen::VectorXd residuals;
    simulation_setup::calculateResiduals< double, double >(
            std::make_shared< ObservationDataset< double, double > >( dataset ), simulators, residuals );

    // With zero-valued simulators, residual computation must return the observed scalar vector itself.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( residuals, expectedObservations, 1.0E-15 );
}

/*!
 * Verifies that legacy observation wrappers remain dataset-backed interfaces.
 *
 * Test outline: constructs SingleObservationSet and ObservationCollection objects
 * through the legacy API, then checks that filtering, restoration, parser-based
 * selection, weight/residual/observation mutation and EstimationInput/
 * CovarianceAnalysisInput references are delegated to the same dataset backend.
 */
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
    const FlattenedObservationData< double, double > singleData = singleDataset->createOrderedFlattenedObservationData( );

    // Converting a legacy single set to a dataset must preserve observations, residuals and weights.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( singleData.getObservationVector( ), angularConversionSet->getObservationsVector( ), 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( singleData.getResidualVector( ), angularConversionSet->getResidualsVector( ), 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( singleData.getWeightVector( ), angularConversionSet->getWeightsVector( ), 1.0E-15 );

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

    // Legacy single-set filtering must reject the selected epoch and expose it through the saved rejected set.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( singleSet->getObservationsVector( ), expectedSingleFilteredVector, 1.0E-15 );
    Eigen::VectorXd expectedSavedFilteredVector( 1 );
    expectedSavedFilteredVector << 11.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            singleSet->getFilteredObservationSet( )->getObservationsVector( ), expectedSavedFilteredVector, 1.0E-15 );

    singleSet->filterObservations( observationFilter( epochs_filtering, std::vector< double >( { 2.0 } ), false ), true );
    Eigen::VectorXd expectedSingleRestoredVector( 3 );
    expectedSingleRestoredVector << 10.0, 11.0, 12.0;

    // Reversing the epoch filter must restore the previously rejected observation.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( singleSet->getObservationsVector( ), expectedSingleRestoredVector, 1.0E-15 );

    singleSet->filterObservations( observationFilter( absolute_value_filtering, 11.0 ), false );
    Eigen::VectorXd expectedSingleDroppedVector( 2 );
    expectedSingleDroppedVector << 10.0, 11.0;

    // Dropping by absolute value without saving rejected rows must update the active legacy view only.
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

    // Collection-to-dataset conversion must use the same ordered flattened data for observations, residuals and weights.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( collectionDataset->createOrderedFlattenedObservationData( ).getObservationVector( ),
                                       observationCollection.getObservationVector( ),
                                       1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( collectionDataset->createOrderedFlattenedObservationData( ).getResidualVector( ),
                                       observationCollection.getConcatenatedResiduals( ),
                                       1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( collectionDataset->createOrderedFlattenedObservationData( ).getWeightVector( ),
                                       observationCollection.getConcatenatedWeights( ),
                                       1.0E-15 );

    Eigen::VectorXd expectedInitialObservationVector( 9 );
    expectedInitialObservationVector << 10.0, 11.0, 20.0, 21.0, 22.0, 23.0, 30.0, 31.0, 32.0;

    // The ordered collection vector must still be observable-type ordered after construction from unsorted input sets.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( observationCollection.getObservationVector( ), expectedInitialObservationVector, 1.0E-15 );

    const std::vector< ObservableType > observableTypes = observationCollection.getObservableTypes( );

    // Collection metadata queries must report the unique observable types represented in the dataset backend.
    BOOST_CHECK_EQUAL( observableTypes.size( ), 3 );
    BOOST_CHECK_EQUAL( std::count( observableTypes.begin( ), observableTypes.end( ), one_way_range ), 1 );
    BOOST_CHECK_EQUAL( std::count( observableTypes.begin( ), observableTypes.end( ), angular_position ), 1 );
    BOOST_CHECK_EQUAL( std::count( observableTypes.begin( ), observableTypes.end( ), position_observable ), 1 );

    const std::vector< LinkEnds > collectionLinkEnds = observationCollection.getLinkEnds( );

    // Link-end metadata must be de-duplicated across sets sharing the same link geometry.
    BOOST_CHECK_EQUAL( collectionLinkEnds.size( ), 2 );
    BOOST_CHECK_EQUAL( std::count( collectionLinkEnds.begin( ), collectionLinkEnds.end( ), station1LinkDefinition.linkEnds_ ), 1 );
    BOOST_CHECK_EQUAL( std::count( collectionLinkEnds.begin( ), collectionLinkEnds.end( ), station2LinkDefinition.linkEnds_ ), 1 );

    const std::vector< std::string > bodiesInLinkEnds = observationCollection.getBodiesInLinkEnds( );

    // Body-name metadata must include each body appearing anywhere in the collection link ends.
    BOOST_CHECK_EQUAL( bodiesInLinkEnds.size( ), 2 );
    BOOST_CHECK_EQUAL( std::count( bodiesInLinkEnds.begin( ), bodiesInLinkEnds.end( ), "Earth" ), 1 );
    BOOST_CHECK_EQUAL( std::count( bodiesInLinkEnds.begin( ), bodiesInLinkEnds.end( ), "Vehicle" ), 1 );

    const std::vector< std::string > referencePointsInLinkEnds = observationCollection.getReferencePointsInLinkEnds( );

    // Reference-point metadata must include each station represented in the link ends.
    BOOST_CHECK_EQUAL( referencePointsInLinkEnds.size( ), 2 );
    BOOST_CHECK_EQUAL( std::count( referencePointsInLinkEnds.begin( ), referencePointsInLinkEnds.end( ), "Station1" ), 1 );
    BOOST_CHECK_EQUAL( std::count( referencePointsInLinkEnds.begin( ), referencePointsInLinkEnds.end( ), "Station2" ), 1 );

    const std::vector< std::pair< double, double > > timeBoundsList = observationCollection.getTimeBoundsList( );

    // Per-set time bounds must be preserved for scalar, angular and position observations.
    BOOST_CHECK_EQUAL( timeBoundsList.size( ), 3 );
    BOOST_CHECK_EQUAL( std::count( timeBoundsList.begin( ), timeBoundsList.end( ), std::make_pair( 1.0, 2.0 ) ), 1 );
    BOOST_CHECK_EQUAL( std::count( timeBoundsList.begin( ), timeBoundsList.end( ), std::make_pair( 3.0, 4.0 ) ), 1 );
    BOOST_CHECK_EQUAL( std::count( timeBoundsList.begin( ), timeBoundsList.end( ), std::make_pair( 5.0, 5.0 ) ), 1 );

    std::vector< std::shared_ptr< SingleObservationSet< double, double > > > rangeObservationSets =
            observationCollection.getSingleObservationSets( observationParser( one_way_range ) );

    // Observable-type parsing must select the same range set exposed by the dataset backend.
    BOOST_REQUIRE_EQUAL( rangeObservationSets.size( ), 1 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            rangeObservationSets.at( 0 )->getObservationsVector( ), rangeSet->getObservationsVector( ), 1.0E-15 );

    std::vector< std::shared_ptr< SingleObservationSet< double, double > > > station1ObservationSets =
            observationCollection.getSingleObservationSets( observationParser( "Station1", true ) );

    // Body/reference-point parsing must find both Station1 sets regardless of observable type.
    BOOST_CHECK_EQUAL( station1ObservationSets.size( ), 2 );

    std::vector< std::shared_ptr< SingleObservationSet< double, double > > > station2ObservationSets =
            observationCollection.getSingleObservationSets( observationParser( station2LinkDefinition.linkEnds_ ) );

    // Exact link-end parsing must isolate the Station2 position set.
    BOOST_REQUIRE_EQUAL( station2ObservationSets.size( ), 1 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            station2ObservationSets.at( 0 )->getObservationsVector( ), positionSet->getObservationsVector( ), 1.0E-15 );

    std::vector< std::shared_ptr< ObservationCollectionParser > > station1RangeParserList;
    station1RangeParserList.push_back( observationParser( one_way_range ) );
    station1RangeParserList.push_back( observationParser( "Station1", true ) );
    std::vector< std::shared_ptr< SingleObservationSet< double, double > > > station1RangeObservationSets =
            observationCollection.getSingleObservationSets( observationParser( station1RangeParserList, true ) );

    // Combined parsers must apply intersection semantics when requested.
    BOOST_REQUIRE_EQUAL( station1RangeObservationSets.size( ), 1 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            station1RangeObservationSets.at( 0 )->getObservationsVector( ), rangeSet->getObservationsVector( ), 1.0E-15 );

    std::vector< std::shared_ptr< SingleObservationSet< double, double > > > notStation1RangeObservationSets =
            observationCollection.getSingleObservationSets( observationParser( station1RangeParserList, true, true ) );

    // Negating a multi-type parser must apply to the combined result, not be ignored by the legacy adapter.
    BOOST_CHECK_EQUAL( notStation1RangeObservationSets.size( ), 2 );

    const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > rangeObservations =
            observationCollection.getObservations( observationParser( one_way_range ) );

    // Parsed observation extraction must return the selected set in legacy vector form.
    BOOST_REQUIRE_EQUAL( rangeObservations.size( ), 1 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( rangeObservations.at( 0 ), rangeSet->getObservationsVector( ), 1.0E-15 );

    const std::vector< std::vector< double > > rangeObservationTimes =
            observationCollection.getObservationTimes( observationParser( one_way_range ) );

    // Parsed time extraction must preserve the event times of the selected set.
    BOOST_REQUIRE_EQUAL( rangeObservationTimes.size( ), 1 );
    const std::vector< double > expectedRangeObservationTimes = rangeSet->getObservationTimes( );
    BOOST_CHECK_EQUAL_COLLECTIONS( rangeObservationTimes.at( 0 ).begin( ),
                                   rangeObservationTimes.at( 0 ).end( ),
                                   expectedRangeObservationTimes.begin( ),
                                   expectedRangeObservationTimes.end( ) );

    const std::map< ObservableType, std::pair< int, int > > observationTypeStartSize =
            observationCollection.getObservationTypeStartAndSize( );

    // Observable-type start/size bookkeeping must match the concatenated scalar-vector layout.
    BOOST_CHECK_EQUAL( observationTypeStartSize.at( one_way_range ).first, 0 );
    BOOST_CHECK_EQUAL( observationTypeStartSize.at( one_way_range ).second, 2 );
    BOOST_CHECK_EQUAL( observationTypeStartSize.at( angular_position ).first, 2 );
    BOOST_CHECK_EQUAL( observationTypeStartSize.at( angular_position ).second, 4 );
    BOOST_CHECK_EQUAL( observationTypeStartSize.at( position_observable ).first, 6 );
    BOOST_CHECK_EQUAL( observationTypeStartSize.at( position_observable ).second, 3 );

    const std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > > observationSetStartSize =
            observationCollection.getObservationSetStartAndSize( );

    // Per-link start/size bookkeeping must resolve each legacy set to its scalar range.
    BOOST_CHECK_EQUAL( observationSetStartSize.at( one_way_range ).at( station1LinkDefinition.linkEnds_ ).at( 0 ).first, 0 );
    BOOST_CHECK_EQUAL( observationSetStartSize.at( one_way_range ).at( station1LinkDefinition.linkEnds_ ).at( 0 ).second, 2 );
    BOOST_CHECK_EQUAL( observationSetStartSize.at( angular_position ).at( station1LinkDefinition.linkEnds_ ).at( 0 ).first, 2 );
    BOOST_CHECK_EQUAL( observationSetStartSize.at( angular_position ).at( station1LinkDefinition.linkEnds_ ).at( 0 ).second, 4 );
    BOOST_CHECK_EQUAL( observationSetStartSize.at( position_observable ).at( station2LinkDefinition.linkEnds_ ).at( 0 ).first, 6 );
    BOOST_CHECK_EQUAL( observationSetStartSize.at( position_observable ).at( station2LinkDefinition.linkEnds_ ).at( 0 ).second, 3 );

    const std::vector< std::pair< int, int > > concatenatedStartSize = observationCollection.getConcatenatedObservationSetStartAndSize( );

    // Concatenated set ranges must expose the same scalar spans without the observable/link grouping.
    BOOST_REQUIRE_EQUAL( concatenatedStartSize.size( ), 3 );
    BOOST_CHECK_EQUAL( concatenatedStartSize.at( 0 ).first, 0 );
    BOOST_CHECK_EQUAL( concatenatedStartSize.at( 0 ).second, 2 );
    BOOST_CHECK_EQUAL( concatenatedStartSize.at( 1 ).first, 2 );
    BOOST_CHECK_EQUAL( concatenatedStartSize.at( 1 ).second, 4 );
    BOOST_CHECK_EQUAL( concatenatedStartSize.at( 2 ).first, 6 );
    BOOST_CHECK_EQUAL( concatenatedStartSize.at( 2 ).second, 3 );

    const std::vector< double > expectedEventTimes = { 1.0, 2.0, 3.0, 4.0, 5.0 };
    const std::vector< double > legacyEventTimes = observationCollection.getConcatenatedObservationTimes( );

    // Legacy event-time vectors must list one time per observation rather than one per scalar component.
    BOOST_CHECK_EQUAL_COLLECTIONS(
            legacyEventTimes.begin( ), legacyEventTimes.end( ), expectedEventTimes.begin( ), expectedEventTimes.end( ) );

    observationCollection.setConstantWeight( 9.0, observationParser( one_way_range ) );
    Eigen::VectorXd expectedWeightsAfterRangeConstant( 9 );
    expectedWeightsAfterRangeConstant << 9.0, 9.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0;

    // Weight mutation through the legacy collection must update both the legacy vector and backing dataset.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( observationCollection.getConcatenatedWeights( ), expectedWeightsAfterRangeConstant, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            observationCollection.getObservationDataset( )->createOrderedFlattenedObservationData( ).getWeightVector( ),
            expectedWeightsAfterRangeConstant,
            1.0E-15 );

    Eigen::VectorXd replacementObservationVector( 9 );
    replacementObservationVector << 100.0, 101.0, 200.0, 201.0, 202.0, 203.0, 300.0, 301.0, 302.0;
    observationCollection.setObservations( replacementObservationVector );

    // Replacing all observations through the legacy API must update both value accessors.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( observationCollection.getObservationVector( ), replacementObservationVector, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( observationCollection.getObservationVectorReference( ), replacementObservationVector, 1.0E-15 );

    const std::vector< double > expectedScalarTimes = { 1.0, 2.0, 3.0, 3.0, 4.0, 4.0, 5.0, 5.0, 5.0 };
    const std::vector< double > legacyScalarTimes = observationCollection.getConcatenatedTimeVector( );

    // Legacy scalar-time vectors must repeat times for angular and position components.
    BOOST_CHECK_EQUAL_COLLECTIONS(
            legacyScalarTimes.begin( ), legacyScalarTimes.end( ), expectedScalarTimes.begin( ), expectedScalarTimes.end( ) );

    observationCollection.filterObservations(
            observationFilter( epochs_filtering, std::vector< double >( { 2.0 } ) ), observationParser( one_way_range ), true );
    Eigen::VectorXd expectedCollectionFilteredVector( 8 );
    expectedCollectionFilteredVector << 100.0, 200.0, 201.0, 202.0, 203.0, 300.0, 301.0, 302.0;

    // Filtering a parsed subset must remove only the matching range epoch from the active collection.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( observationCollection.getObservationVector( ), expectedCollectionFilteredVector, 1.0E-15 );

    observationCollection.filterObservations(
            observationFilter( epochs_filtering, std::vector< double >( { 2.0 } ), false ), observationParser( one_way_range ), true );

    // Reversing the parsed filter must restore the complete legacy observation vector.
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

    // Estimation and covariance inputs created from a legacy collection must share the collection dataset.
    BOOST_CHECK( covarianceInputDataset == estimationCollection->getObservationDataset( ) );
    BOOST_CHECK( estimationInputDataset == estimationCollection->getObservationDataset( ) );

    Eigen::VectorXd residuals( 3 );
    residuals << 0.1, 0.2, 0.3;
    estimationCollection->setResiduals( residuals );

    // Residual updates through the legacy collection must be visible through both input objects.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            covarianceInputDataset->createOrderedFlattenedObservationData( ).getResidualVector( ), residuals, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            estimationInputDataset->createOrderedFlattenedObservationData( ).getResidualVector( ), residuals, 1.0E-15 );

    estimationCollection->filterObservations(
            observationFilter( epochs_filtering, std::vector< double >( { 2.0 } ) ), observationParser( one_way_range ), true );

    // Filtering must update the shared dataset in place rather than replacing the input-object references.
    BOOST_CHECK( covarianceInputDataset == estimationCollection->getObservationDataset( ) );
    BOOST_CHECK( estimationInputDataset == estimationCollection->getObservationDataset( ) );
    BOOST_CHECK_EQUAL( covarianceInputDataset->getTotalScalarSize( ), 2 );
    BOOST_CHECK_EQUAL( estimationInputDataset->getTotalScalarSize( ), 2 );
    BOOST_CHECK_EQUAL( estimationCollection->getTotalObservableSize( ), 2 );

    const std::vector< double > expectedTimes = { 1.0, 3.0 };
    const std::vector< double > covarianceInputTimes = covarianceInputDataset->createOrderedFlattenedObservationData( ).getTimes( );
    const std::vector< double > estimationInputTimes = estimationInputDataset->createOrderedFlattenedObservationData( ).getTimes( );
    const std::vector< double > collectionTimes = estimationCollection->getConcatenatedTimeVector( );

    // After filtering, the remaining times must agree across covariance input, estimation input and collection views.
    BOOST_CHECK_EQUAL_COLLECTIONS(
            covarianceInputTimes.begin( ), covarianceInputTimes.end( ), expectedTimes.begin( ), expectedTimes.end( ) );
    BOOST_CHECK_EQUAL_COLLECTIONS(
            estimationInputTimes.begin( ), estimationInputTimes.end( ), expectedTimes.begin( ), expectedTimes.end( ) );
    BOOST_CHECK_EQUAL_COLLECTIONS( collectionTimes.begin( ), collectionTimes.end( ), expectedTimes.begin( ), expectedTimes.end( ) );

    estimationCollection->filterObservations(
            observationFilter( epochs_filtering, std::vector< double >( { 2.0 } ), false ), observationParser( one_way_range ), true );

    // Restoring observations must again keep the input objects tied to the same full dataset.
    BOOST_CHECK( covarianceInputDataset == estimationCollection->getObservationDataset( ) );
    BOOST_CHECK( estimationInputDataset == estimationCollection->getObservationDataset( ) );
    BOOST_CHECK_EQUAL( covarianceInputDataset->getTotalScalarSize( ), 3 );
    BOOST_CHECK_EQUAL( estimationInputDataset->getTotalScalarSize( ), 3 );
    const std::vector< double > expectedRestoredTimes = { 1.0, 2.0, 3.0 };
    const std::vector< double > restoredCovarianceInputTimes = covarianceInputDataset->createOrderedFlattenedObservationData( ).getTimes( );
    const std::vector< double > restoredEstimationInputTimes = estimationInputDataset->createOrderedFlattenedObservationData( ).getTimes( );

    // Restored input flattened data must recover the original event-time sequence.
    BOOST_CHECK_EQUAL_COLLECTIONS( restoredCovarianceInputTimes.begin( ),
                                   restoredCovarianceInputTimes.end( ),
                                   expectedRestoredTimes.begin( ),
                                   expectedRestoredTimes.end( ) );
    BOOST_CHECK_EQUAL_COLLECTIONS( restoredEstimationInputTimes.begin( ),
                                   restoredEstimationInputTimes.end( ),
                                   expectedRestoredTimes.begin( ),
                                   expectedRestoredTimes.end( ) );

    covarianceInputDataset->addObservationSet(
            one_way_range, station1LinkDefinition, { ( Eigen::VectorXd( 1 ) << 13.0 ).finished( ) }, { 4.0 }, receiver );

    // Direct dataset mutation through a shared input pointer must refresh the legacy collection caches on next access.
    BOOST_CHECK_EQUAL( estimationCollection->getTotalObservableSize( ), covarianceInputDataset->getTotalScalarSize( ) );
    BOOST_CHECK_EQUAL( estimationCollection->getObservationsSets( ).at( one_way_range ).at( station1LinkDefinition.linkEnds_ ).size( ), 2 );
    const std::vector< double > externallyMutatedTimes = estimationCollection->getConcatenatedTimeVector( );
    const std::vector< double > expectedExternallyMutatedTimes = { 1.0, 2.0, 3.0, 4.0 };

    // The refreshed collection flattened data must include rows added outside the legacy wrapper.
    BOOST_CHECK_EQUAL_COLLECTIONS( externallyMutatedTimes.begin( ),
                                   externallyMutatedTimes.end( ),
                                   expectedExternallyMutatedTimes.begin( ),
                                   expectedExternallyMutatedTimes.end( ) );
}

BOOST_AUTO_TEST_CASE( test_legacy_weight_setters_delegate_to_dataset_backend )
{
    const LinkDefinition linkDefinition = createOneWayLinkDefinition( "Station1" );

    const auto createManualObservationSet =
            []( const ObservableType observableType,
                const LinkDefinition& linkDefinition,
                const std::vector< Eigen::VectorXd >& observations,
                const std::vector< double >& times ) -> std::shared_ptr< SingleObservationSet< double, double > > {
        return std::make_shared< SingleObservationSet< double, double > >( observableType, linkDefinition, observations, times, receiver );
    };

    std::shared_ptr< SingleObservationSet< double, double > > rangeSet = createManualObservationSet(
            one_way_range,
            linkDefinition,
            { Eigen::Vector1d::Constant( 10.0 ), Eigen::Vector1d::Constant( 20.0 ), Eigen::Vector1d::Constant( 30.0 ) },
            { 1.0, 2.0, 3.0 } );
    std::shared_ptr< SingleObservationSet< double, double > > angularSet = createManualObservationSet(
            angular_position,
            linkDefinition,
            { ( Eigen::Vector2d( ) << 40.0, 41.0 ).finished( ), ( Eigen::Vector2d( ) << 50.0, 51.0 ).finished( ) },
            { 4.0, 5.0 } );

    ObservationCollection< double, double > observationCollection(
            std::vector< std::shared_ptr< SingleObservationSet< double, double > > >( { angularSet, rangeSet } ) );
    std::shared_ptr< ObservationDataset< double, double > > dataset = observationCollection.getObservationDataset( );

    const auto checkWeights = [ & ]( const Eigen::VectorXd& expectedWeights ) {
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( observationCollection.getConcatenatedWeights( ), expectedWeights, 1.0E-15 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( dataset->createOrderedFlattenedObservationData( ).getWeightVector( ), expectedWeights, 1.0E-15 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( dataset->createOrderedFlattenedObservationData( ).getSparseWeightMatrix( ).toDense( ),
                                           expectedWeights.asDiagonal( ),
                                           1.0E-15 );
    };

    std::map< std::shared_ptr< ObservationCollectionParser >, double > constantWeights;
    constantWeights[ observationParser( one_way_range ) ] = 4.0;
    constantWeights[ observationParser( angular_position ) ] = 9.0;
    observationCollection.setConstantWeightPerObservable( constantWeights );

    Eigen::VectorXd expectedConstantWeights( 7 );
    expectedConstantWeights << 4.0, 4.0, 4.0, 9.0, 9.0, 9.0, 9.0;

    // Parser-map constant weights must update the legacy vector and the backing dataset in the same ordered layout.
    checkWeights( expectedConstantWeights );

    Eigen::VectorXd rangeTabulatedWeights( 3 );
    rangeTabulatedWeights << 1.1, 1.2, 1.3;
    Eigen::VectorXd angularTabulatedWeights( 4 );
    angularTabulatedWeights << 2.1, 2.2, 2.3, 2.4;

    std::map< std::shared_ptr< ObservationCollectionParser >, Eigen::VectorXd > tabulatedWeights;
    tabulatedWeights[ observationParser( one_way_range ) ] = rangeTabulatedWeights;
    tabulatedWeights[ observationParser( angular_position ) ] = angularTabulatedWeights;
    observationCollection.setTabulatedWeights( tabulatedWeights );

    Eigen::VectorXd expectedTabulatedWeights( 7 );
    expectedTabulatedWeights << 1.1, 1.2, 1.3, 2.1, 2.2, 2.3, 2.4;

    // Parser-map tabulated weights must remain diagonal and agree with the dataset flattened sparse matrix.
    checkWeights( expectedTabulatedWeights );

    const ObservationDatasetViewer< double, double > staleViewer =
            dataset->createViewer( ObservationSelectionCondition< double, double >::all( ) );
    observationCollection.removeSingleObservationSets( observationParser( one_way_range ) );

    // Removing sets through the legacy collection must structurally mutate the shared dataset and invalidate old viewers.
    BOOST_CHECK_THROW( staleViewer.getNumberOfObservations( ), std::runtime_error );
    BOOST_CHECK_EQUAL( dataset->getNumberOfObservationSets( ), 1 );
    BOOST_CHECK_EQUAL( dataset->getObservationSetMetadata( 0 ).observableType_, angular_position );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( observationCollection.getConcatenatedWeights( ), angularTabulatedWeights, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            dataset->createOrderedFlattenedObservationData( ).getWeightVector( ), angularTabulatedWeights, 1.0E-15 );
}

/*!
 * Verifies the legacy ObservationCollection observation-set splitters.
 *
 * Test outline: builds small legacy collections and applies every splitter type
 * still exposed through ObservationCollection: explicit time tags, maximum time
 * span, maximum number of observations and maximum adjacent time interval. The
 * assertions check both the number of generated sets and the row/time ranges in
 * each set, including the legacy minimum-observation rule that drops too-small
 * split fragments.
 */
BOOST_AUTO_TEST_CASE( test_legacy_observation_collection_splitters )
{
    const LinkDefinition station1LinkDefinition = createOneWayLinkDefinition( "Station1" );
    const LinkDefinition station2LinkDefinition = createOneWayLinkDefinition( "Station2" );

    const auto createScalarSet = []( const LinkDefinition& linkDefinition,
                                     const std::vector< double >& times ) -> std::shared_ptr< SingleObservationSet< double, double > > {
        std::vector< Eigen::VectorXd > observations;
        for( std::size_t i = 0; i < times.size( ); ++i )
        {
            // The scalar value records the original row number, making it easy
            // to verify that split sets preserve row order and row association.
            observations.push_back( ( Eigen::VectorXd( 1 ) << static_cast< double >( i ) ).finished( ) );
        }
        return std::make_shared< SingleObservationSet< double, double > >( one_way_range, linkDefinition, observations, times, receiver );
    };

    const auto getRangeSetsForStation1 = []( ObservationCollection< double, double >& collection )
            -> std::vector< std::shared_ptr< SingleObservationSet< double, double > > > {
        return collection.getSingleObservationSets( observationParser( "Station1", true ) );
    };

    const auto checkTimes = []( const std::vector< double >& actualTimes, const std::vector< double >& expectedTimes ) {
        // Each split set must expose exactly the expected ordered event times.
        BOOST_CHECK_EQUAL_COLLECTIONS( actualTimes.begin( ), actualTimes.end( ), expectedTimes.begin( ), expectedTimes.end( ) );
    };

    {
        ObservationCollection< double, double > collection( std::vector< std::shared_ptr< SingleObservationSet< double, double > > >(
                { createScalarSet( station1LinkDefinition, { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 } ) } ) );

        collection.splitObservationSets( observationSetSplitter( time_tags_splitter, std::vector< double >( { 0.0, 3.0, 6.0 } ) ),
                                         observationParser( "Station1", true ) );

        const std::vector< std::shared_ptr< SingleObservationSet< double, double > > > splitSets = getRangeSetsForStation1( collection );

        // Time-tag splitting should create one set up to the split tag and one after it.
        BOOST_REQUIRE_EQUAL( splitSets.size( ), 2 );
        BOOST_CHECK_EQUAL( splitSets.at( 0 )->getNumberOfObservables( ), 4 );
        BOOST_CHECK_EQUAL( splitSets.at( 1 )->getNumberOfObservables( ), 3 );
        checkTimes( splitSets.at( 0 )->getObservationTimes( ), { 0.0, 1.0, 2.0, 3.0 } );
        checkTimes( splitSets.at( 1 )->getObservationTimes( ), { 4.0, 5.0, 6.0 } );
    }

    {
        ObservationCollection< double, double > collection( std::vector< std::shared_ptr< SingleObservationSet< double, double > > >(
                { createScalarSet( station1LinkDefinition, { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 } ) } ) );

        collection.splitObservationSets( observationSetSplitter( time_span_splitter, 3.0, 3 ), observationParser( "Station1", true ) );

        const std::vector< std::shared_ptr< SingleObservationSet< double, double > > > splitSets = getRangeSetsForStation1( collection );

        // The final two-row fragment is below minNumberObservations and should be dropped by the legacy splitter.
        BOOST_REQUIRE_EQUAL( splitSets.size( ), 2 );
        checkTimes( splitSets.at( 0 )->getObservationTimes( ), { 0.0, 1.0, 2.0, 3.0 } );
        checkTimes( splitSets.at( 1 )->getObservationTimes( ), { 4.0, 5.0, 6.0, 7.0 } );
        for( const std::shared_ptr< SingleObservationSet< double, double > >& splitSet : splitSets )
        {
            // Every retained split set must satisfy the requested maximum span.
            BOOST_CHECK( splitSet->getObservationTimes( ).back( ) - splitSet->getObservationTimes( ).front( ) <= 3.0 );
        }
    }

    {
        ObservationCollection< double, double > collection( std::vector< std::shared_ptr< SingleObservationSet< double, double > > >(
                { createScalarSet( station1LinkDefinition, { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 } ) } ) );

        collection.splitObservationSets( observationSetSplitter( nb_observations_splitter, 4, 4 ), observationParser( "Station1", true ) );

        const std::vector< std::shared_ptr< SingleObservationSet< double, double > > > splitSets = getRangeSetsForStation1( collection );

        // Splitting by maximum count should keep two complete four-row sets and drop the two-row tail.
        BOOST_REQUIRE_EQUAL( splitSets.size( ), 2 );
        BOOST_CHECK_EQUAL( splitSets.at( 0 )->getNumberOfObservables( ), 4 );
        BOOST_CHECK_EQUAL( splitSets.at( 1 )->getNumberOfObservables( ), 4 );
        checkTimes( splitSets.at( 0 )->getObservationTimes( ), { 0.0, 1.0, 2.0, 3.0 } );
        checkTimes( splitSets.at( 1 )->getObservationTimes( ), { 4.0, 5.0, 6.0, 7.0 } );
    }

    {
        ObservationCollection< double, double > collection( std::vector< std::shared_ptr< SingleObservationSet< double, double > > >(
                { createScalarSet( station1LinkDefinition, { 0.0, 1.0, 2.0, 8.0, 9.0 } ),
                  createScalarSet( station2LinkDefinition, { 0.0, 100.0 } ) } ) );

        collection.splitObservationSets( observationSetSplitter( time_interval_splitter, 2.0 ), observationParser( "Station1", true ) );

        const std::vector< std::shared_ptr< SingleObservationSet< double, double > > > station1SplitSets =
                getRangeSetsForStation1( collection );
        const std::vector< std::shared_ptr< SingleObservationSet< double, double > > > station2Sets =
                collection.getSingleObservationSets( observationParser( "Station2", true ) );

        // A large adjacent time gap in Station1 should create two sets.
        BOOST_REQUIRE_EQUAL( station1SplitSets.size( ), 2 );
        checkTimes( station1SplitSets.at( 0 )->getObservationTimes( ), { 0.0, 1.0, 2.0 } );
        checkTimes( station1SplitSets.at( 1 )->getObservationTimes( ), { 8.0, 9.0 } );

        // The parser passed to splitObservationSets must keep unrelated Station2 data untouched.
        BOOST_REQUIRE_EQUAL( station2Sets.size( ), 1 );
        checkTimes( station2Sets.at( 0 )->getObservationTimes( ), { 0.0, 100.0 } );
    }
}

BOOST_AUTO_TEST_CASE( test_legacy_collection_preserves_single_set_sharing_and_dependent_variable_clearing )
{
    const LinkDefinition linkDefinition = createOneWayLinkDefinition( "Station1" );
    std::shared_ptr< SingleObservationSet< double, double > > sharedSet = std::make_shared< SingleObservationSet< double, double > >(
            one_way_range,
            linkDefinition,
            std::vector< Eigen::VectorXd >( { Eigen::Vector1d::Constant( 1.0 ), Eigen::Vector1d::Constant( 2.0 ) } ),
            std::vector< double >( { 1.0, 2.0 } ),
            receiver );

    ObservationCollection< double, double > collection(
            std::vector< std::shared_ptr< SingleObservationSet< double, double > > >( { sharedSet } ) );
    const ObservationCollection< double, double >::SortedObservationSets observationSets = collection.getObservationsSets( );
    BOOST_CHECK( observationSets.at( one_way_range ).at( linkDefinition.linkEnds_ ).at( 0 ) == sharedSet );

    std::shared_ptr< SingleObservationSet< double, double > > sortedConstructorSet =
            std::make_shared< SingleObservationSet< double, double > >(
                    one_way_range,
                    linkDefinition,
                    std::vector< Eigen::VectorXd >( { Eigen::Vector1d::Constant( 5.0 ), Eigen::Vector1d::Constant( 6.0 ) } ),
                    std::vector< double >( { 3.0, 4.0 } ),
                    receiver );
    ObservationCollection< double, double >::SortedObservationSets sortedInputSets;
    sortedInputSets[ one_way_range ][ linkDefinition.linkEnds_ ].push_back( sortedConstructorSet );
    ObservationCollection< double, double > sortedConstructorCollection( sortedInputSets );
    BOOST_CHECK( sortedConstructorCollection.getObservationsSets( ).at( one_way_range ).at( linkDefinition.linkEnds_ ).at( 0 ) ==
                 sortedConstructorSet );

    sharedSet->setObservations( std::vector< Eigen::VectorXd >( { Eigen::Vector1d::Constant( 3.0 ), Eigen::Vector1d::Constant( 4.0 ) } ) );
    Eigen::VectorXd expectedMutatedObservations( 2 );
    expectedMutatedObservations << 3.0, 4.0;

    // Mutating the original SingleObservationSet after collection construction must still update the collection view.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( collection.getObservationVector( ), expectedMutatedObservations, 1.0E-15 );

    std::vector< Eigen::VectorXd > dependentVariables(
            { ( Eigen::Vector2d( ) << 0.1, 0.2 ).finished( ), ( Eigen::Vector2d( ) << 0.3, 0.4 ).finished( ) } );
    sharedSet->setObservationsDependentVariables( dependentVariables );
    BOOST_CHECK_EQUAL( sharedSet->getObservationsDependentVariables( ).size( ), 2 );

    std::vector< Eigen::VectorXd > emptyDependentVariables;
    sharedSet->setObservationsDependentVariables( emptyDependentVariables );

    // The legacy empty-vector setter is the documented clearing path and must not throw.
    BOOST_CHECK( sharedSet->getObservationsDependentVariables( ).empty( ) );
    BOOST_CHECK( collection.getObservationsSets( ).at( one_way_range ).at( linkDefinition.linkEnds_ ).at( 0 ) == sharedSet );
}

BOOST_AUTO_TEST_CASE( test_weighted_design_matrix_output_uses_sparse_weights )
{
    Eigen::MatrixXd normalizedDesignMatrix( 2, 2 );
    normalizedDesignMatrix << 1.0, 2.0, 3.0, 4.0;
    Eigen::Matrix2d denseWeights;
    denseWeights << 4.0, 1.0, 1.0, 3.0;
    const Eigen::SparseMatrix< double > sparseWeights = denseWeights.sparseView( );
    const Eigen::Vector2d diagonalWeights = denseWeights.diagonal( );
    const Eigen::Vector2d normalizationTerms = Eigen::Vector2d::Ones( );

    simulation_setup::CovarianceAnalysisOutput< double, double > covarianceOutput( normalizedDesignMatrix,
                                                                                   diagonalWeights,
                                                                                   normalizationTerms,
                                                                                   Eigen::Matrix2d::Identity( ),
                                                                                   Eigen::MatrixXd::Zero( 0, 0 ),
                                                                                   Eigen::VectorXd::Zero( 0 ),
                                                                                   Eigen::MatrixXd::Zero( 0, 0 ),
                                                                                   Eigen::MatrixXd::Zero( 0, 0 ),
                                                                                   false,
                                                                                   sparseWeights );

    const Eigen::MatrixXd lowerWeightSquareRoot = Eigen::LLT< Eigen::MatrixXd >( denseWeights ).matrixL( );
    const Eigen::MatrixXd expectedWeightedDesignMatrix = lowerWeightSquareRoot.transpose( ) * normalizedDesignMatrix;

    BOOST_CHECK( covarianceOutput.hasFullWeightMatrix( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( covarianceOutput.getWeightsMatrix( ).toDense( ), denseWeights, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( covarianceOutput.getNormalizedWeightedDesignMatrix( ), expectedWeightedDesignMatrix, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( covarianceOutput.getUnnormalizedWeightedDesignMatrix( ), expectedWeightedDesignMatrix, 1.0E-15 );
}

BOOST_AUTO_TEST_CASE( test_large_sparse_weighted_design_matrix_uses_sparse_cholesky )
{
    const int numberOfObservations = 5000;
    const int firstCorrelatedRow = 123;
    const int secondCorrelatedRow = 4321;
    const double diagonalWeight = 4.0;
    const double offDiagonalWeight = 1.0;

    Eigen::MatrixXd normalizedDesignMatrix = Eigen::MatrixXd::Zero( numberOfObservations, 2 );
    normalizedDesignMatrix.row( firstCorrelatedRow ) << 1.0, 2.0;
    normalizedDesignMatrix.row( secondCorrelatedRow ) << -3.0, 0.5;
    normalizedDesignMatrix.row( 2500 ) << 0.25, -1.5;

    std::vector< Eigen::Triplet< double > > weightTriplets;
    weightTriplets.reserve( numberOfObservations + 2 );
    for( int i = 0; i < numberOfObservations; ++i )
    {
        weightTriplets.emplace_back( i, i, diagonalWeight );
    }
    weightTriplets.emplace_back( firstCorrelatedRow, secondCorrelatedRow, offDiagonalWeight );
    weightTriplets.emplace_back( secondCorrelatedRow, firstCorrelatedRow, offDiagonalWeight );

    Eigen::SparseMatrix< double > sparseWeights( numberOfObservations, numberOfObservations );
    sparseWeights.setFromTriplets( weightTriplets.begin( ), weightTriplets.end( ) );
    const Eigen::VectorXd diagonalWeights = Eigen::VectorXd::Constant( numberOfObservations, diagonalWeight );
    const Eigen::Vector2d normalizationTerms = Eigen::Vector2d::Ones( );

    simulation_setup::CovarianceAnalysisOutput< double, double > covarianceOutput( normalizedDesignMatrix,
                                                                                   diagonalWeights,
                                                                                   normalizationTerms,
                                                                                   Eigen::Matrix2d::Identity( ),
                                                                                   Eigen::MatrixXd::Zero( 0, 0 ),
                                                                                   Eigen::VectorXd::Zero( 0 ),
                                                                                   Eigen::MatrixXd::Zero( 0, 0 ),
                                                                                   Eigen::MatrixXd::Zero( 0, 0 ),
                                                                                   false,
                                                                                   sparseWeights );

    BOOST_CHECK( covarianceOutput.hasFullWeightMatrix( ) );
    BOOST_CHECK( !covarianceOutput.isSparseWeightCholeskyFactorCurrent( ) );

    const Eigen::MatrixXd weightedDesignMatrix = covarianceOutput.getNormalizedWeightedDesignMatrix( );

    BOOST_CHECK( covarianceOutput.isSparseWeightCholeskyFactorCurrent( ) );
    BOOST_CHECK_EQUAL( covarianceOutput.getWeightsMatrix( ).nonZeros( ), numberOfObservations + 2 );

    Eigen::MatrixXd expectedWeightedRows( 3, 2 );
    expectedWeightedRows.row( 0 ) = std::sqrt( diagonalWeight ) * normalizedDesignMatrix.row( firstCorrelatedRow ) +
            offDiagonalWeight / std::sqrt( diagonalWeight ) * normalizedDesignMatrix.row( secondCorrelatedRow );
    expectedWeightedRows.row( 1 ) = std::sqrt( diagonalWeight - offDiagonalWeight * offDiagonalWeight / diagonalWeight ) *
            normalizedDesignMatrix.row( secondCorrelatedRow );
    expectedWeightedRows.row( 2 ) = std::sqrt( diagonalWeight ) * normalizedDesignMatrix.row( 2500 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( weightedDesignMatrix.row( firstCorrelatedRow ), expectedWeightedRows.row( 0 ), 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( weightedDesignMatrix.row( secondCorrelatedRow ), expectedWeightedRows.row( 1 ), 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( weightedDesignMatrix.row( 2500 ), expectedWeightedRows.row( 2 ), 1.0E-15 );
}

/*!
 * Verifies empty-set behavior and invalid input validation.
 *
 * Test outline: creates an empty observation set and checks that all counters,
 * flattened data and time bounds are well-defined. It then exercises inconsistent
 * observation/time/weight dimensions and invalid set mutation requests to
 * ensure they fail with runtime errors.
 */
BOOST_AUTO_TEST_CASE( test_dataset_empty_sets_and_invalid_inputs )
{
    const LinkDefinition linkDefinition = createOneWayLinkDefinition( "Station1" );
    ObservationDataset< double, double > dataset;

    const int emptySetId = dataset.addObservationSet(
            one_way_range, linkDefinition, std::vector< Eigen::VectorXd >( ), std::vector< double >( ), receiver );

    // Empty sets must remain valid metadata containers with zero-length flattened data and undefined time bounds.
    BOOST_CHECK_EQUAL( dataset.getNumberOfObservationSets( ), 1 );
    BOOST_CHECK_EQUAL( dataset.getNumberOfObservationsForSet( emptySetId ), 0 );
    BOOST_CHECK_EQUAL( dataset.getTotalScalarSizeForSet( emptySetId ), 0 );
    BOOST_CHECK_EQUAL( dataset.createEstimationFlattenedObservationData( ).getObservationVector( ).size( ), 0 );
    BOOST_CHECK( std::isnan( dataset.getTimeBoundsForSet( emptySetId ).first ) );
    BOOST_CHECK( std::isnan( dataset.getTimeBoundsForSet( emptySetId ).second ) );

    dataset.setLinkEndReferencePoint( "Earth", "RenamedStation", transmitter, ObservationSelectionCondition< double, double >::all( ) );
    BOOST_CHECK_EQUAL( dataset.getLinkDefinition( dataset.getObservationSetMetadata( emptySetId ).linkDefinitionId_ )
                               .linkEnds_.at( transmitter )
                               .getReferencePointName( ),
                       "Station1" );

    // Inconsistent dimensions or mutating a nonexistent row in an empty set must be rejected.
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
    BOOST_CHECK_THROW(
            dataset.addObservationSet( angular_position, linkDefinition, { Eigen::Vector1d::Constant( 1.0 ) }, { 1.0 }, receiver ),
            std::runtime_error );

    const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > elevationSettings =
            simulation_setup::elevationAngleDependentVariable( receiver, LinkEndId( "Vehicle", "" ) );
    const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping > bookkeeping =
            std::make_shared< simulation_setup::ObservationDependentVariableBookkeeping >( one_way_range, linkDefinition );
    bookkeeping->addDependentVariable( elevationSettings );
    BOOST_CHECK_THROW( dataset.addObservationSet( one_way_range,
                                                  linkDefinition,
                                                  { Eigen::Vector1d::Constant( 1.0 ), Eigen::Vector1d::Constant( 2.0 ) },
                                                  { 1.0, 2.0 },
                                                  receiver,
                                                  { Eigen::Vector1d::Constant( 3.0 ), Eigen::Vector2d::Constant( 4.0 ) },
                                                  bookkeeping ),
                       std::runtime_error );

    const unsigned int layoutOnlySetId = dataset.addObservationSet( one_way_range,
                                                                    linkDefinition,
                                                                    { Eigen::Vector1d::Constant( 1.0 ) },
                                                                    { 1.0 },
                                                                    receiver,
                                                                    std::vector< Eigen::VectorXd >( ),
                                                                    bookkeeping );
    BOOST_CHECK_THROW( dataset.getSingleDependentVariableForSet( layoutOnlySetId, elevationSettings ), std::runtime_error );
    BOOST_CHECK( dataset.getAllCompatibleDependentVariablesForSet( layoutOnlySetId, elevationSettings ).empty( ) );
    BOOST_CHECK_THROW( dataset.setDependentVariablesForSet( layoutOnlySetId, { Eigen::Vector2d::Ones( ) } ), std::runtime_error );
    BOOST_CHECK_THROW( dataset.addObservationsToSet(
                               layoutOnlySetId, { Eigen::Vector1d::Constant( 2.0 ) }, { 2.0 }, { Eigen::Vector1d::Constant( 3.0 ) } ),
                       std::runtime_error );
    BOOST_CHECK_THROW( dataset.setObservationVectorForSet( emptySetId, ( Eigen::VectorXd( 1 ) << 1.0 ).finished( ) ), std::runtime_error );
    BOOST_CHECK_THROW( dataset.removeObservationsFromSet( emptySetId, std::vector< unsigned int >( { 0 } ) ), std::runtime_error );
}

/*!
 * Verifies sorting, duplicate removal and observation transfer edge cases.
 *
 * Test outline: adds intentionally unsorted and duplicate observations, removes
 * duplicates, validates sorted row order, copies and moves selected observations
 * between sets and confirms invalid transfer requests are rejected.
 */
BOOST_AUTO_TEST_CASE( test_dataset_duplicate_selection_and_move_edge_cases )
{
    const LinkDefinition linkDefinition = createOneWayLinkDefinition( "Station1" );

    ObservationDataset< double, double > dataset;
    const unsigned int sourceSetId = dataset.addObservationSet( one_way_range,
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

    // Duplicate removal must keep one row at a duplicated epoch and sort the retained rows by time.
    BOOST_CHECK_EQUAL( dataset.getNumberOfObservationsForSet( sourceSetId ), 3 );
    const std::vector< double > sortedUniqueTimes = dataset.getObservationTimesForSet( sourceSetId );
    const std::vector< double > expectedSortedUniqueTimes = { 1.0, 3.0, 4.0 };
    BOOST_CHECK_EQUAL_COLLECTIONS(
            sortedUniqueTimes.begin( ), sortedUniqueTimes.end( ), expectedSortedUniqueTimes.begin( ), expectedSortedUniqueTimes.end( ) );

    checkIds( dataset.getObservationIdsMatchingCondition( ObservationSelectionCondition< double, double >::residualAbsoluteValueGreaterThan(
                      ( Eigen::VectorXd( 1 ) << 0.25 ).finished( ) ) ),
              { 1 } );

    // Value- and residual-based conditions must operate on the sorted, duplicate-free rows.
    checkIds( dataset.getObservationIdsMatchingCondition(
                      ObservationSelectionCondition< double, double >::observationAbsoluteValueGreaterThan(
                              ( Eigen::VectorXd( 1 ) << 12.0 ).finished( ) ) ),
              { 1 } );

    ObservationDataset< double, double > targetDataset;
    const unsigned int targetSetId = targetDataset.addObservationSet(
            one_way_range, linkDefinition, std::vector< Eigen::VectorXd >( ), std::vector< double >( ), receiver );

    dataset.moveObservationsToSet( sourceSetId, targetDataset, targetSetId, std::vector< unsigned int >( { 0, 2 } ), false );

    // Copying rows to another dataset must leave the source set untouched.
    BOOST_CHECK_EQUAL( dataset.getNumberOfObservationsForSet( sourceSetId ), 3 );
    BOOST_CHECK_EQUAL( targetDataset.getNumberOfObservationsForSet( targetSetId ), 2 );

    const ObservationDatasetViewer< double, double > staleMoveViewer =
            dataset.createViewer( ObservationSelectionCondition< double, double >::all( ) );
    dataset.moveObservationsToSet( sourceSetId, targetDataset, targetSetId, std::vector< unsigned int >( { 1 } ), true );

    // Moving rows must shrink the source set while appending the moved observation to the target set.
    BOOST_CHECK_THROW( staleMoveViewer.getNumberOfObservations( ), std::runtime_error );
    BOOST_CHECK_EQUAL( dataset.getNumberOfObservationsForSet( sourceSetId ), 2 );
    BOOST_CHECK_EQUAL( targetDataset.getNumberOfObservationsForSet( targetSetId ), 3 );

    const Eigen::VectorXd targetObservationVector = targetDataset.getObservationVectorForSet( targetSetId );
    Eigen::VectorXd expectedTargetObservationVector( 3 );
    expectedTargetObservationVector << 11.0, 13.0, 10.0;

    // Copied and moved observations must retain their values and append order in the target dataset.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( targetObservationVector, expectedTargetObservationVector, 1.0E-15 );

    ObservationDataset< double, double > duplicateDataset;
    const unsigned int duplicateSetId = duplicateDataset.addObservationSet(
            one_way_range,
            linkDefinition,
            { Eigen::Vector1d::Constant( 1.0 ), Eigen::Vector1d::Constant( 2.0 ), Eigen::Vector1d::Constant( 3.0 ) },
            { 1.0, 2.0, 1.0 },
            receiver );
    const ObservationDatasetViewer< double, double > staleDuplicateViewer =
            duplicateDataset.createViewer( ObservationSelectionCondition< double, double >::all( ) );
    std::ostringstream duplicateWarningStream;
    std::streambuf* originalWarningBuffer = std::cerr.rdbuf( duplicateWarningStream.rdbuf( ) );
    duplicateDataset.eraseDuplicateObservationsFromSet( duplicateSetId, false );
    std::cerr.rdbuf( originalWarningBuffer );

    // Direct duplicate erasure must remove duplicate epochs even when they are not adjacent, optionally suppress warnings, and invalidate
    // old viewers.
    BOOST_CHECK_EQUAL( duplicateDataset.getNumberOfObservationsForSet( duplicateSetId ), 2 );
    BOOST_CHECK( duplicateWarningStream.str( ).empty( ) );
    BOOST_CHECK_THROW( staleDuplicateViewer.getNumberOfObservations( ), std::runtime_error );

    ObservationDataset< double, double > warningDuplicateDataset;
    const unsigned int warningDuplicateSetId = warningDuplicateDataset.addObservationSet(
            one_way_range,
            linkDefinition,
            { Eigen::Vector1d::Constant( 1.0 ), Eigen::Vector1d::Constant( 2.0 ), Eigen::Vector1d::Constant( 3.0 ) },
            { 1.0, 1.0, 2.0 },
            receiver );
    duplicateWarningStream.str( "" );
    duplicateWarningStream.clear( );
    originalWarningBuffer = std::cerr.rdbuf( duplicateWarningStream.rdbuf( ) );
    warningDuplicateDataset.eraseDuplicateObservationsFromSet( warningDuplicateSetId, true );
    std::cerr.rdbuf( originalWarningBuffer );

    BOOST_CHECK( duplicateWarningStream.str( ).find( "Detected and removed 1" ) != std::string::npos );
    BOOST_CHECK( duplicateWarningStream.str( ).find( "duplicate observations" ) != std::string::npos );
}

/*!
 * Verifies row-level ObservationSelectionCondition selection logic.
 *
 * Test outline: creates observations with link metadata, values, residuals, weights,
 * dependent variables and active/rejected status. It checks that each condition
 * selects exactly the intended observation ids and that combined conditions can
 * be used to build consistent viewers.
 */
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

    // Set-id conditions must be backend conditions, not Python-only convenience filters.
    const ObservationSelectionCondition< double, double > secondSet = ObservationSelectionCondition< double, double >::setId( 1 );
    BOOST_CHECK( secondSet.getConditionType( ) == ObservationSelectionConditionType::set_id );
    BOOST_CHECK_EQUAL( secondSet.getConditionTypeString( ), "set_id" );
    BOOST_CHECK_EQUAL( secondSet.getSetId( ), 1 );
    checkIds( dataset.getObservationIdsMatchingCondition( secondSet ), { 3, 4 } );

    // Link-based conditions must match rows by complete link definition, link-end type and individual link end.
    checkIds( dataset.getObservationIdsMatchingCondition(
                      ObservationSelectionCondition< double, double >::linkDefinition( station2LinkDefinition ) ),
              { 3, 4 } );
    checkIds( dataset.getObservationIdsMatchingCondition( ObservationSelectionCondition< double, double >::linkEndType( transmitter ) ),
              { 0, 1, 2, 3, 4 } );
    checkIds( dataset.getObservationIdsMatchingCondition(
                      ObservationSelectionCondition< double, double >::linkEnd( transmitter, LinkEndId( "Earth", "Station1" ) ) ),
              { 0, 1, 2 } );

    const ObservationSelectionCondition< double, double > rangeOrLastAngular =
            ObservationSelectionCondition< double, double >::observableType( one_way_range ) ||
            ObservationSelectionCondition< double, double >::timeBounds( 4.5, 5.5 );

    // Query-tree metadata must expose the logical operator and its two typed child conditions.
    BOOST_CHECK( rangeOrLastAngular.getConditionType( ) == ObservationSelectionConditionType::or_condition );
    BOOST_CHECK_EQUAL( rangeOrLastAngular.getChildConditions( ).size( ), 2 );
    BOOST_CHECK( rangeOrLastAngular.getChildConditions( ).at( 0 ).getConditionType( ) ==
                 ObservationSelectionConditionType::observable_type );
    BOOST_CHECK( rangeOrLastAngular.getChildConditions( ).at( 1 ).getConditionType( ) == ObservationSelectionConditionType::time_bounds );
    BOOST_CHECK_EQUAL( rangeOrLastAngular.getChildConditions( ).at( 0 ).getObservableType( ), one_way_range );
    BOOST_CHECK_CLOSE_FRACTION( rangeOrLastAngular.getChildConditions( ).at( 1 ).getTimeBounds( ).first, 4.5, 1.0E-15 );

    // One-sided time comparison conditions must expose their comparison value and select expected rows.
    const ObservationSelectionCondition< double, double > afterOrAtFour =
            ObservationSelectionCondition< double, double >::timeGreaterEqual( 4.0 );
    BOOST_CHECK( afterOrAtFour.getConditionType( ) == ObservationSelectionConditionType::time_greater_equal );
    BOOST_CHECK_CLOSE_FRACTION( afterOrAtFour.getTimeValue( ), 4.0, 1.0E-15 );
    checkIds( dataset.getObservationIdsMatchingCondition( afterOrAtFour ), { 3, 4 } );
    checkIds( dataset.getObservationIdsMatchingCondition( ObservationSelectionCondition< double, double >::timeGreaterThan( 4.0 ) ),
              { 4 } );
    checkIds( dataset.getObservationIdsMatchingCondition( ObservationSelectionCondition< double, double >::timeLessEqual( 2.0 ) ),
              { 0, 1 } );
    checkIds( dataset.getObservationIdsMatchingCondition( ObservationSelectionCondition< double, double >::timeLessThan( 2.0 ) ), { 0 } );

    // Boolean condition composition must support unions and complements of row selections.
    checkIds( dataset.getObservationIdsMatchingCondition( rangeOrLastAngular ), { 0, 1, 2, 4 } );

    const ObservationSelectionCondition< double, double > outsideMiddleTimes =
            !ObservationSelectionCondition< double, double >::timeBounds( 2.0, 4.0 );

    // Negated conditions must expose their child while evaluating the complement selection.
    BOOST_CHECK( outsideMiddleTimes.getConditionType( ) == ObservationSelectionConditionType::not_condition );
    BOOST_CHECK_EQUAL( outsideMiddleTimes.getChildConditions( ).size( ), 1 );
    checkIds( dataset.getObservationIdsMatchingCondition( outsideMiddleTimes ), { 0, 4 } );

    // Value, residual and dependent-variable conditions must compare component-wise vectors of matching size.
    checkIds(
            dataset.getObservationIdsMatchingCondition( ObservationSelectionCondition< double, double >::observableType( one_way_range ) &&
                                                        ObservationSelectionCondition< double, double >::residualAbsoluteValueGreaterThan(
                                                                ( Eigen::VectorXd( 1 ) << 0.5 ).finished( ) ) ),
            { 1 } );

    checkIds( dataset.getObservationIdsMatchingCondition(
                      ObservationSelectionCondition< double, double >::observableType( angular_position ) &&
                      ObservationSelectionCondition< double, double >::observationAbsoluteValueGreaterThan(
                              ( Eigen::Vector2d( ) << 2.5, 3.5 ).finished( ) ) ),
              { 4 } );

    checkIds( dataset.getObservationIdsMatchingCondition( ObservationSelectionCondition< double, double >::dependentVariableGreaterThan(
                      elevationSettings, ( Eigen::VectorXd( 1 ) << 0.3 ).finished( ) ) ),
              { 1 } );

    // Scalar limits broadcast across vector-valued observations and select a row when any component exceeds the limit.
    checkIds( dataset.getObservationIdsMatchingCondition(
                      ObservationSelectionCondition< double, double >::observableType( angular_position ) &&
                      ObservationSelectionCondition< double, double >::residualAbsoluteValueGreaterThan( 0.6 ) ),
              { 4 } );
    checkIds( dataset.getObservationIdsMatchingCondition(
                      ObservationSelectionCondition< double, double >::observableType( angular_position ) &&
                      ObservationSelectionCondition< double, double >::observationAbsoluteValueGreaterThan( 3.5 ) ),
              { 4 } );
    checkIds( dataset.getObservationIdsMatchingCondition(
                      ObservationSelectionCondition< double, double >::dependentVariableGreaterThan( elevationSettings, 0.3 ) ),
              { 1 } );

    std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping > ambiguousDependentVariableBookkeeping =
            std::make_shared< simulation_setup::ObservationDependentVariableBookkeeping >( one_way_range, station1LinkDefinition );
    ambiguousDependentVariableBookkeeping->addDependentVariable( elevationSettings );
    ambiguousDependentVariableBookkeeping->addDependentVariable( elevationSettings );

    ObservationDataset< double, double > ambiguousDependentVariableDataset;
    ambiguousDependentVariableDataset.addObservationSet(
            one_way_range,
            station1LinkDefinition,
            { ( Eigen::VectorXd( 1 ) << 10.0 ).finished( ), ( Eigen::VectorXd( 1 ) << 20.0 ).finished( ) },
            { 1.0, 2.0 },
            receiver,
            { ( Eigen::Vector2d( ) << 0.1, 10.0 ).finished( ), ( Eigen::Vector2d( ) << 0.4, 20.0 ).finished( ) },
            ambiguousDependentVariableBookkeeping );

    BOOST_CHECK_THROW(
            ambiguousDependentVariableDataset.getObservationIdsMatchingCondition(
                    ObservationSelectionCondition< double, double >::dependentVariableGreaterThan( elevationSettings, 0.3, false ) ),
            std::runtime_error );
    checkIds( ambiguousDependentVariableDataset.getObservationIdsMatchingCondition(
                      ObservationSelectionCondition< double, double >::dependentVariableGreaterThan( elevationSettings, 0.3, true ) ),
              { 1 } );

    ObservationDataset< double, double > noDependentVariableBookkeepingDataset;
    noDependentVariableBookkeepingDataset.addObservationSet(
            one_way_range,
            station1LinkDefinition,
            { ( Eigen::VectorXd( 1 ) << 10.0 ).finished( ), ( Eigen::VectorXd( 1 ) << 20.0 ).finished( ) },
            { 1.0, 2.0 },
            receiver );
    checkIds( noDependentVariableBookkeepingDataset.getObservationIdsMatchingCondition(
                      ObservationSelectionCondition< double, double >::dependentVariableGreaterThan( elevationSettings, 0.3 ) ),
              {} );

    // Conditions with incompatible vector sizes must fail before selecting rows.
    BOOST_CHECK_THROW(
            dataset.getObservationIdsMatchingCondition( ObservationSelectionCondition< double, double >::residualAbsoluteValueGreaterThan(
                    ( Eigen::Vector2d( ) << 0.1, 0.2 ).finished( ) ) ),
            std::runtime_error );
    BOOST_CHECK_THROW(
            dataset.getObservationIdsMatchingCondition( ObservationSelectionCondition< double, double >::dependentVariableGreaterThan(
                    elevationSettings, ( Eigen::Vector2d( ) << 0.1, 0.2 ).finished( ) ) ),
            std::runtime_error );

    dataset.rejectObservations( ObservationSelectionCondition< double, double >::observableType( one_way_range ) &&
                                        ObservationSelectionCondition< double, double >::residualAbsoluteValueGreaterThan(
                                                ( Eigen::VectorXd( 1 ) << 0.5 ).finished( ) ),
                                "large residual" );

    // Rejection and restoration conditions must toggle row status without changing row ids.
    checkIds( dataset.getObservationIdsMatchingCondition( ObservationSelectionCondition< double, double >::rejected( ) ), { 1 } );
    checkIds( dataset.getObservationIdsMatchingCondition( ObservationSelectionCondition< double, double >::active( ) ), { 0, 2, 3, 4 } );
    dataset.restoreObservations( ObservationSelectionCondition< double, double >::rejected( ) );
    checkIds( dataset.getObservationIdsMatchingCondition( ObservationSelectionCondition< double, double >::active( ) ), { 0, 1, 2, 3, 4 } );
}

/*!
 * Verifies that viewer ordered flattening follows dataset ordered-output order.
 *
 * Test outline: inserts angular-position rows before one-way range rows, so
 * dataset row order differs from ordered flattened-data order. A viewer
 * selecting all rows must keep dataset order for estimation flattening but use
 * observable/link/set ordered output for ordered flattening.
 */
BOOST_AUTO_TEST_CASE( test_dataset_viewer_ordered_flattening_reorders_selected_rows )
{
    ObservationDataset< double, double > dataset;
    const LinkDefinition station1LinkDefinition = createOneWayLinkDefinition( "Station1" );
    const LinkDefinition station2LinkDefinition = createOneWayLinkDefinition( "Station2" );

    dataset.addObservationSet( angular_position,
                               station2LinkDefinition,
                               { ( Eigen::Vector2d( ) << 20.0, 21.0 ).finished( ), ( Eigen::Vector2d( ) << 22.0, 23.0 ).finished( ) },
                               { 3.0, 4.0 },
                               receiver );
    dataset.addObservationSet( one_way_range,
                               station1LinkDefinition,
                               { ( Eigen::VectorXd( 1 ) << 10.0 ).finished( ), ( Eigen::VectorXd( 1 ) << 11.0 ).finished( ) },
                               { 1.0, 2.0 },
                               receiver );

    const ObservationDatasetViewer< double, double > viewer =
            dataset.createViewer( ObservationSelectionCondition< double, double >::all( ) );

    // The viewer selection itself is in dataset row order: angular rows first, then range rows.
    checkIds( viewer.getObservationIds( ), { 0, 1, 2, 3 } );

    // Estimation flattening intentionally preserves the viewer's selected row order.
    checkIds( viewer.createEstimationFlattenedObservationData( ).getObservationIds( ), { 0, 0, 1, 1, 2, 3 } );

    // Ordered flattening must instead follow legacy ordered output: one-way range before angular position.
    checkIds( viewer.createOrderedFlattenedObservationData( ).getObservationIds( ), { 2, 3, 0, 0, 1, 1 } );

    std::map< ObservableType, std::shared_ptr< ObservationSimulatorBase< double, double > > > simulators;
    simulators[ one_way_range ] = std::make_shared< ZeroObservationSimulator >( one_way_range, 1 );
    simulators[ angular_position ] = std::make_shared< ZeroObservationSimulator >( angular_position, 2 );
    Eigen::VectorXd legacyResiduals;
    simulation_setup::calculateResiduals< double, double >(
            createObservationCollection( std::make_shared< ObservationDataset< double, double > >( dataset ) ),
            simulators,
            legacyResiduals );
    Eigen::VectorXd expectedLegacyResiduals( 6 );
    expectedLegacyResiduals << 10.0, 11.0, 20.0, 21.0, 22.0, 23.0;

    // Legacy collections expose residuals in their observable/link ordered layout, independent of dataset insertion order.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( legacyResiduals, expectedLegacyResiduals, 1.0E-15 );
}

/*!
 * Verifies viewers, rejection/restoration and reduced dataset creation.
 *
 * Test outline: builds a dataset, selects rows with conditions, creates viewers and
 * reduced datasets, rejects/restores observations and checks flattened data sizes.
 * It also confirms that structural mutations invalidate previously created
 * viewers.
 */
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

    const ObservationSelectionCondition< double, double > middleTimes =
            ObservationSelectionCondition< double, double >::observableType( one_way_range ) &&
            ObservationSelectionCondition< double, double >::timeBounds( 1.5, 3.5 );
    const std::vector< unsigned int > selectedObservationIds = dataset.getObservationIdsMatchingCondition( middleTimes );
    const std::vector< unsigned int > expectedSelectedObservationIds = { 1, 2 };

    // A time-window condition must identify the middle observations by their stable row ids.
    BOOST_CHECK_EQUAL_COLLECTIONS( selectedObservationIds.begin( ),
                                   selectedObservationIds.end( ),
                                   expectedSelectedObservationIds.begin( ),
                                   expectedSelectedObservationIds.end( ) );

    const ObservationDatasetViewer< double, double > viewer = dataset.createViewer( middleTimes );

    // A viewer must project only rows that satisfy its condition.
    BOOST_CHECK_EQUAL( viewer.getNumberOfObservations( ), 2 );
    Eigen::VectorXd expectedViewerObservations( 2 );
    expectedViewerObservations << 20.0, 30.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            viewer.createEstimationFlattenedObservationData( ).getObservationVector( ), expectedViewerObservations, 1.0E-15 );

    const std::shared_ptr< ObservationDataset< double, double > > keptDataset = dataset.createNewAndKeep( middleTimes );

    // A reduced keep-dataset must contain only selected rows and preserve their scalar flattened data.
    BOOST_CHECK_EQUAL( keptDataset->getNumberOfObservationSets( ), 1 );
    BOOST_CHECK_EQUAL( keptDataset->getNumberOfObservations( ), 2 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            keptDataset->createEstimationFlattenedObservationData( ).getObservationVector( ), expectedViewerObservations, 1.0E-15 );

    const std::shared_ptr< ObservationDataset< double, double > > droppedDataset = dataset.createNewAndDrop( middleTimes );
    Eigen::VectorXd expectedDroppedObservations( 2 );
    expectedDroppedObservations << 10.0, 40.0;

    // A reduced drop-dataset must contain the complement of the selected rows.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            droppedDataset->createEstimationFlattenedObservationData( ).getObservationVector( ), expectedDroppedObservations, 1.0E-15 );

    dataset.rejectObservations( ObservationSelectionCondition< double, double >::timeBounds( 2.5, 3.5 ), "test rejection" );

    // Rejection must store row status and reason while removing the row from active flattened data.
    BOOST_CHECK( !dataset.getObservationRow( 2 ).isActive_ );
    BOOST_CHECK_EQUAL( dataset.getObservationRow( 2 ).rejectionReason_, "test rejection" );

    Eigen::VectorXd expectedActiveObservations( 3 );
    expectedActiveObservations << 10.0, 20.0, 40.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            dataset.createEstimationFlattenedObservationData( ).getObservationVector( ), expectedActiveObservations, 1.0E-15 );
    BOOST_CHECK_EQUAL( dataset.createEstimationFlattenedObservationData( true ).getObservationVector( ).size( ), 4 );
    BOOST_CHECK_EQUAL( dataset.createViewer( ObservationSelectionCondition< double, double >::rejected( ) ).getNumberOfObservations( ), 1 );

    dataset.restoreObservations( ObservationSelectionCondition< double, double >::rejected( ) );

    // Restoration must reactivate rejected rows and clear their rejection reason.
    BOOST_CHECK( dataset.getObservationRow( 2 ).isActive_ );
    BOOST_CHECK( dataset.getObservationRow( 2 ).rejectionReason_.empty( ) );
    BOOST_CHECK_EQUAL( dataset.createEstimationFlattenedObservationData( ).getObservationVector( ).size( ), 4 );

    dataset.rejectObservations( ObservationSelectionCondition< double, double >::timeBounds( 2.5, 3.5 ), "delete rejection" );
    dataset.removeRejectedObservations( );

    // Deleting rejected rows is explicit and physically removes rows from the dataset.
    Eigen::VectorXd expectedAfterRejectedDeletion( 3 );
    expectedAfterRejectedDeletion << 10.0, 20.0, 40.0;
    BOOST_CHECK_EQUAL( dataset.getNumberOfObservations( ), 3 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            dataset.createEstimationFlattenedObservationData( true ).getObservationVector( ), expectedAfterRejectedDeletion, 1.0E-15 );
    BOOST_CHECK_EQUAL( dataset.createViewer( ObservationSelectionCondition< double, double >::rejected( ) ).getNumberOfObservations( ), 0 );

    dataset.removeObservations( ObservationSelectionCondition< double, double >::timeBounds( 3.5, 4.5 ) );

    // Direct condition-based removal must provide the destructive counterpart to createNewAndDrop.
    Eigen::VectorXd expectedAfterConditionRemoval( 2 );
    expectedAfterConditionRemoval << 10.0, 20.0;
    BOOST_CHECK_EQUAL( dataset.getNumberOfObservations( ), 2 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            dataset.createEstimationFlattenedObservationData( true ).getObservationVector( ), expectedAfterConditionRemoval, 1.0E-15 );

    const ObservationDatasetViewer< double, double > invalidatedViewer =
            dataset.createViewer( ObservationSelectionCondition< double, double >::all( ) );
    dataset.addObservationSet( one_way_range, linkDefinition, std::vector< Eigen::VectorXd >( ), std::vector< double >( ), receiver );

    // Structural dataset mutations must invalidate previously created viewers.
    BOOST_CHECK_THROW( invalidatedViewer.getNumberOfObservations( ), std::runtime_error );

    ObservationDataset< double, double > linkMutationDataset;
    linkMutationDataset.addObservationSet( one_way_range, linkDefinition, { Eigen::Vector1d::Constant( 1.0 ) }, { 1.0 }, receiver );
    const ObservationDatasetViewer< double, double > linkMutationViewer =
            linkMutationDataset.createViewer( ObservationSelectionCondition< double, double >::all( ) );
    linkMutationDataset.setLinkEndReferencePoint( "Earth", "StationX", transmitter );

    // Updating link-end reference points changes set metadata and must invalidate old viewers.
    BOOST_CHECK_THROW( linkMutationViewer.getNumberOfObservations( ), std::runtime_error );
    BOOST_CHECK_EQUAL( linkMutationDataset.getLinkDefinition( linkMutationDataset.getObservationSetMetadata( 0 ).linkDefinitionId_ )
                               .linkEnds_.at( transmitter )
                               .getReferencePointName( ),
                       "StationX" );
}

BOOST_AUTO_TEST_CASE( test_dataset_viewer_lifetime_and_legacy_cache_invalidation )
{
    const LinkDefinition station1LinkDefinition = createOneWayLinkDefinition( "Station1" );
    const LinkDefinition station2LinkDefinition = createOneWayLinkDefinition( "Station2" );

    const ObservationDatasetViewer< double, double > expiredViewer = [ station1LinkDefinition ]( ) {
        ObservationDataset< double, double > localDataset;
        localDataset.addObservationSet( one_way_range, station1LinkDefinition, { Eigen::Vector1d::Constant( 1.0 ) }, { 1.0 }, receiver );
        return localDataset.createViewer( ObservationSelectionCondition< double, double >::all( ) );
    }( );
    BOOST_CHECK_THROW( expiredViewer.getNumberOfObservations( ), std::runtime_error );

    ObservationDataset< double, double > assignedDataset;
    assignedDataset.addObservationSet( one_way_range, station1LinkDefinition, { Eigen::Vector1d::Constant( 1.0 ) }, { 1.0 }, receiver );
    const ObservationDatasetViewer< double, double > preAssignmentViewer =
            assignedDataset.createViewer( ObservationSelectionCondition< double, double >::all( ) );
    ObservationDataset< double, double > replacementDataset;
    replacementDataset.addObservationSet( one_way_range, station1LinkDefinition, { Eigen::Vector1d::Constant( 2.0 ) }, { 2.0 }, receiver );
    assignedDataset = replacementDataset;
    BOOST_CHECK_THROW( preAssignmentViewer.getNumberOfObservations( ), std::runtime_error );

    const std::shared_ptr< ObservationDataset< double, double > > collectionDataset =
            std::make_shared< ObservationDataset< double, double > >( );
    collectionDataset->addObservationSet( one_way_range, station1LinkDefinition, { Eigen::Vector1d::Constant( 1.0 ) }, { 1.0 }, receiver );
    ObservationCollection< double, double > collection( collectionDataset );
    collection.getSingleObservationSets( ).front( )->resetLinkEnds( station2LinkDefinition );

    const std::vector< LinkDefinition > refreshedLinkDefinitions = collection.getLinkDefinitionsForSingleObservable( one_way_range );
    BOOST_REQUIRE_EQUAL( refreshedLinkDefinitions.size( ), 1 );
    BOOST_CHECK( refreshedLinkDefinitions.front( ) == station2LinkDefinition );
}

/*!
 * Verifies compact, block and sparse off-diagonal weight handling.
 *
 * Test outline: covers scalar weights, per-observation blocks, set-level blocks,
 * observation-id selected cross blocks, symmetric transpose insertion,
 * component selection and rejection/restoration of weighted observations. It
 * checks both compact weight vectors and materialized sparse flattened data
 * matrices.
 */
BOOST_AUTO_TEST_CASE( test_dataset_compact_and_matrix_weights )
{
    const LinkDefinition linkDefinition = createOneWayLinkDefinition( "Station1" );

    ObservationDataset< double, double > scalarWeightDataset;
    const int scalarSetId = scalarWeightDataset.addObservationSetWithWeights(
            angular_position,
            linkDefinition,
            { ( Eigen::Vector2d( ) << 10.0, 11.0 ).finished( ), ( Eigen::Vector2d( ) << 20.0, 21.0 ).finished( ) },
            { 1.0, 2.0 },
            receiver,
            ObservationWeightSettings::constantScalar( 5.0 ) );

    // Scalar-weight construction must remain diagonal-only and expose no per-observation matrix blocks.
    BOOST_CHECK_EQUAL( scalarSetId, 0 );
    Eigen::VectorXd expectedScalarWeights( 4 );
    expectedScalarWeights << 5.0, 5.0, 5.0, 5.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( scalarWeightDataset.getWeightVectorForSet( scalarSetId ), expectedScalarWeights, 1.0E-15 );
    const Eigen::MatrixXd expectedScalarWeightMatrix = 5.0 * Eigen::MatrixXd::Identity( 4, 4 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( scalarWeightDataset.createEstimationFlattenedObservationData( ).getSparseWeightMatrix( ).toDense( ),
                                       expectedScalarWeightMatrix,
                                       1.0E-15 );
    BOOST_CHECK_EQUAL( scalarWeightDataset.createEstimationFlattenedObservationData( ).getSparseWeightMatrix( ).nonZeros( ), 4 );
    BOOST_CHECK( scalarWeightDataset.createEstimationFlattenedObservationData( ).isDiagonalWeightOnly( ) );
    BOOST_CHECK( !scalarWeightDataset.hasWeightMatrixForObservation( scalarWeightDataset.getObservationIdsForSet( scalarSetId ).at( 0 ) ) );

    Eigen::Matrix2d observationWeightBlock;
    observationWeightBlock << 2.0, 0.5, 0.5, 3.0;
    ObservationDataset< double, double > blockWeightDataset;
    const int blockSetId = blockWeightDataset.addObservationSetWithWeights(
            angular_position,
            linkDefinition,
            { ( Eigen::Vector2d( ) << 10.0, 11.0 ).finished( ), ( Eigen::Vector2d( ) << 20.0, 21.0 ).finished( ) },
            { 1.0, 2.0 },
            receiver,
            ObservationWeightSettings::constantBlock( observationWeightBlock ) );
    Eigen::MatrixXd expectedObservationBlockMatrix = Eigen::MatrixXd::Zero( 4, 4 );
    expectedObservationBlockMatrix.block( 0, 0, 2, 2 ) = observationWeightBlock;
    expectedObservationBlockMatrix.block( 2, 2, 2, 2 ) = observationWeightBlock;

    // Per-observation weight blocks must be repeated for each vector-valued observation and force sparse weighting.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( blockWeightDataset.createEstimationFlattenedObservationData( ).getSparseWeightMatrix( ).toDense( ),
                                       expectedObservationBlockMatrix,
                                       1.0E-15 );
    BOOST_CHECK_EQUAL( blockWeightDataset.createEstimationFlattenedObservationData( ).getSparseWeightMatrix( ).nonZeros( ), 8 );
    BOOST_CHECK( blockWeightDataset.createEstimationFlattenedObservationData( ).hasOffDiagonalWeights( ) );
    BOOST_CHECK( blockWeightDataset.hasWeightMatrixForObservation( blockWeightDataset.getObservationIdsForSet( blockSetId ).at( 0 ) ) );
    Eigen::VectorXd expectedBlockWeightDiagonal( 4 );
    expectedBlockWeightDiagonal << 2.0, 3.0, 2.0, 3.0;

    // Compact weight-vector access must return the diagonal of repeated per-observation blocks.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( blockWeightDataset.getWeightVectorForSet( blockSetId ), expectedBlockWeightDiagonal, 1.0E-15 );

    ObservationDataset< double, double > setConstantMatrixDataset;
    const int constantMatrixSetId = setConstantMatrixDataset.addObservationSet(
            angular_position,
            linkDefinition,
            { ( Eigen::Vector2d( ) << 10.0, 11.0 ).finished( ), ( Eigen::Vector2d( ) << 20.0, 21.0 ).finished( ) },
            { 1.0, 2.0 },
            receiver );
    setConstantMatrixDataset.setConstantSingleObservationMatrixWeightForSet( constantMatrixSetId, observationWeightBlock );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            setConstantMatrixDataset.createEstimationFlattenedObservationData( ).getSparseWeightMatrix( ).toDense( ),
            expectedObservationBlockMatrix,
            1.0E-15 );
    BOOST_CHECK_THROW(
            setConstantMatrixDataset.setConstantSingleObservationMatrixWeightForSet( constantMatrixSetId, Eigen::Matrix3d::Identity( ) ),
            std::runtime_error );

    ObservationDataset< double, double > conditionWeightDataset;
    const int conditionRangeSetId =
            conditionWeightDataset.addObservationSet( one_way_range,
                                                      linkDefinition,
                                                      { Eigen::Vector1d::Constant( 10.0 ), Eigen::Vector1d::Constant( 20.0 ) },
                                                      { 1.0, 2.0 },
                                                      receiver );
    const int conditionAngularSetId = conditionWeightDataset.addObservationSet(
            angular_position,
            linkDefinition,
            { ( Eigen::Vector2d( ) << 30.0, 31.0 ).finished( ), ( Eigen::Vector2d( ) << 40.0, 41.0 ).finished( ) },
            { 3.0, 4.0 },
            receiver );
    conditionWeightDataset.setConstantSingleObservationScalarWeight(
            ObservationSelectionCondition< double, double >::observableType( one_way_range ), 7.0 );
    conditionWeightDataset.setConstantSingleObservationDiagonalWeight(
            ObservationSelectionCondition< double, double >::timeBounds( 3.5, 4.5 ), ( Eigen::Vector2d( ) << 8.0, 9.0 ).finished( ) );
    conditionWeightDataset.setConstantSingleObservationMatrixWeight(
            ObservationSelectionCondition< double, double >::timeBounds( 2.5, 3.5 ), observationWeightBlock );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            conditionWeightDataset.getWeightVectorForSet( conditionRangeSetId ), ( 7.0 * Eigen::Vector2d::Ones( ) ), 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( conditionWeightDataset.getWeightMatrixForObservation(
                                               conditionWeightDataset.getObservationIdsForSet( conditionAngularSetId ).at( 0 ) ),
                                       observationWeightBlock,
                                       1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            conditionWeightDataset.getWeightValue( conditionWeightDataset.getObservationIdsForSet( conditionAngularSetId ).at( 1 ) ),
            ( Eigen::Vector2d( ) << 8.0, 9.0 ).finished( ),
            1.0E-15 );
    BOOST_CHECK_THROW( conditionWeightDataset.setConstantSingleObservationDiagonalWeight(
                               ObservationSelectionCondition< double, double >::all( ), ( Eigen::Vector2d( ) << 1.0, 2.0 ).finished( ) ),
                       std::runtime_error );

    Eigen::MatrixXd setWeightBlock = Eigen::MatrixXd::Zero( 4, 4 );
    setWeightBlock << 1.0, 0.1, 0.2, 0.3, 0.1, 2.0, 0.4, 0.5, 0.2, 0.4, 3.0, 0.6, 0.3, 0.5, 0.6, 4.0;
    ObservationDataset< double, double > setBlockWeightDataset;
    const int setBlockSetId = setBlockWeightDataset.addObservationSetWithWeights(
            angular_position,
            linkDefinition,
            { ( Eigen::Vector2d( ) << 10.0, 11.0 ).finished( ), ( Eigen::Vector2d( ) << 20.0, 21.0 ).finished( ) },
            { 1.0, 2.0 },
            receiver,
            ObservationWeightSettings::setBlock( setWeightBlock ) );

    // A set-level weight block must be stored as one matrix spanning the full scalar range of the set.
    BOOST_CHECK( setBlockWeightDataset.hasWeightMatrixForSet( setBlockSetId ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( setBlockWeightDataset.getWeightMatrixForSet( setBlockSetId ), setWeightBlock, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            setBlockWeightDataset.createEstimationFlattenedObservationData( ).getSparseWeightMatrix( ).toDense( ),
            setWeightBlock,
            1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            setBlockWeightDataset.createEstimationFlattenedObservationData( ).getWeightVector( ), setWeightBlock.diagonal( ), 1.0E-15 );

    ObservationDataset< double, double > extraBlockDataset;
    extraBlockDataset.addObservationSetWithWeights(
            one_way_range,
            linkDefinition,
            { Eigen::Vector1d::Constant( 10.0 ), Eigen::Vector1d::Constant( 20.0 ), Eigen::Vector1d::Constant( 30.0 ) },
            { 1.0, 2.0, 3.0 },
            receiver,
            ObservationWeightSettings::constantScalar( 1.0 ) );
    const std::vector< unsigned int > extraBlockObservationIds = extraBlockDataset.getObservationIdsForSet( 0 );
    const Eigen::Matrix2d extraWeightBlock = ( Eigen::Matrix2d( ) << 4.0, 0.25, 0.25, 5.0 ).finished( );
    extraBlockDataset.setWeightBlock( { extraBlockObservationIds.at( 0 ), extraBlockObservationIds.at( 2 ) },
                                      { extraBlockObservationIds.at( 0 ), extraBlockObservationIds.at( 2 ) },
                                      extraWeightBlock,
                                      {},
                                      {} );

    // Arbitrary observation-id blocks must be stored separately and materialized into the sparse flattened data.
    BOOST_CHECK( extraBlockDataset.hasExtraWeightBlocks( ) );
    BOOST_CHECK_EQUAL( extraBlockDataset.getExtraWeightBlocks( ).size( ), 1 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( extraBlockDataset.getExtraWeightBlocks( ).at( 0 ).weightBlock_, extraWeightBlock, 1.0E-15 );
    Eigen::MatrixXd expectedExtraBlockWeightMatrix = Eigen::MatrixXd::Identity( 3, 3 );
    expectedExtraBlockWeightMatrix( 0, 0 ) = 4.0;
    expectedExtraBlockWeightMatrix( 2, 2 ) = 5.0;
    expectedExtraBlockWeightMatrix( 0, 2 ) = 0.25;
    expectedExtraBlockWeightMatrix( 2, 0 ) = 0.25;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( extraBlockDataset.createEstimationFlattenedObservationData( ).getSparseWeightMatrix( ).toDense( ),
                                       expectedExtraBlockWeightMatrix,
                                       1.0E-15 );

    ObservationDataset< double, double > symmetricComponentBlockDataset;
    symmetricComponentBlockDataset.addObservationSetWithWeights( angular_position,
                                                                 linkDefinition,
                                                                 { ( Eigen::Vector2d( ) << 10.0, 11.0 ).finished( ),
                                                                   ( Eigen::Vector2d( ) << 20.0, 21.0 ).finished( ),
                                                                   ( Eigen::Vector2d( ) << 30.0, 31.0 ).finished( ) },
                                                                 { 1.0, 2.0, 3.0 },
                                                                 receiver,
                                                                 ObservationWeightSettings::constantScalar( 1.0 ) );
    const std::vector< unsigned int > symmetricObservationIds = symmetricComponentBlockDataset.getObservationIdsForSet( 0 );
    const Eigen::MatrixXd crossComponentBlock = ( Eigen::Matrix< double, 2, 1 >( ) << 0.7, 0.8 ).finished( );
    symmetricComponentBlockDataset.setWeightBlock( { symmetricObservationIds.at( 0 ), symmetricObservationIds.at( 1 ) },
                                                   { symmetricObservationIds.at( 2 ) },
                                                   crossComponentBlock,
                                                   { 0 },
                                                   { 1 } );

    // Component-selected blocks with symmetric insertion must create both the original and transposed blocks.
    BOOST_CHECK_EQUAL( symmetricComponentBlockDataset.getExtraWeightBlocks( ).size( ), 2 );
    Eigen::MatrixXd expectedSymmetricComponentWeightMatrix = Eigen::MatrixXd::Identity( 6, 6 );
    expectedSymmetricComponentWeightMatrix( 0, 5 ) = 0.7;
    expectedSymmetricComponentWeightMatrix( 2, 5 ) = 0.8;
    expectedSymmetricComponentWeightMatrix( 5, 0 ) = 0.7;
    expectedSymmetricComponentWeightMatrix( 5, 2 ) = 0.8;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            symmetricComponentBlockDataset.createEstimationFlattenedObservationData( ).getSparseWeightMatrix( ).toDense( ),
            expectedSymmetricComponentWeightMatrix,
            1.0E-15 );

    // Component indices outside the observable size must be rejected before a malformed block can be stored.
    BOOST_CHECK_THROW(
            symmetricComponentBlockDataset.setWeightBlock(
                    { symmetricObservationIds.at( 0 ) }, { symmetricObservationIds.at( 1 ) }, Eigen::MatrixXd::Ones( 1, 1 ), { 3 }, { 0 } ),
            std::runtime_error );
    const std::size_t extraWeightBlockCountBeforeInvalidSymmetricBlock = symmetricComponentBlockDataset.getExtraWeightBlocks( ).size( );
    BOOST_CHECK_THROW( symmetricComponentBlockDataset.setWeightBlock( { symmetricObservationIds.at( 0 ), symmetricObservationIds.at( 1 ) },
                                                                      { symmetricObservationIds.at( 0 ), symmetricObservationIds.at( 1 ) },
                                                                      ( Eigen::Matrix2d( ) << 1.0, 2.0, 3.0, 4.0 ).finished( ),
                                                                      {},
                                                                      {} ),
                       std::runtime_error );
    BOOST_CHECK_EQUAL( symmetricComponentBlockDataset.getExtraWeightBlocks( ).size( ), extraWeightBlockCountBeforeInvalidSymmetricBlock );

    ObservationDataset< double, double > invalidAddWeightDataset;
    BOOST_CHECK_THROW( invalidAddWeightDataset.addObservationSetWithWeights(
                               angular_position,
                               linkDefinition,
                               { ( Eigen::Vector2d( ) << 10.0, 11.0 ).finished( ) },
                               { 1.0 },
                               receiver,
                               ObservationWeightSettings::constantBlock( Eigen::Matrix3d::Identity( ) ) ),
                       std::runtime_error );
    BOOST_CHECK_EQUAL( invalidAddWeightDataset.getNumberOfObservationSets( ), 0 );

    BOOST_CHECK_THROW(
            invalidAddWeightDataset.addObservationSetWithWeights(
                    angular_position,
                    linkDefinition,
                    { ( Eigen::Vector2d( ) << 10.0, 11.0 ).finished( ) },
                    { 1.0 },
                    receiver,
                    ObservationWeightSettings::blockPerObservation( std::vector< Eigen::MatrixXd >( { Eigen::Matrix3d::Identity( ) } ) ) ),
            std::runtime_error );
    BOOST_CHECK_EQUAL( invalidAddWeightDataset.getNumberOfObservationSets( ), 0 );

    BOOST_CHECK_THROW(
            invalidAddWeightDataset.addObservationSetWithWeights( angular_position,
                                                                  linkDefinition,
                                                                  { ( Eigen::Vector2d( ) << 10.0, 11.0 ).finished( ) },
                                                                  { 1.0 },
                                                                  receiver,
                                                                  ObservationWeightSettings::setBlock( Eigen::Matrix3d::Identity( ) ) ),
            std::runtime_error );
    BOOST_CHECK_EQUAL( invalidAddWeightDataset.getNumberOfObservationSets( ), 0 );

    ObservationDataset< double, double > rejectedSetBlockDataset;
    Eigen::MatrixXd fullSetWeightBlock = Eigen::MatrixXd::Zero( 6, 6 );
    for( int row = 0; row < fullSetWeightBlock.rows( ); ++row )
    {
        for( int column = 0; column < fullSetWeightBlock.cols( ); ++column )
        {
            fullSetWeightBlock( row, column ) =
                    row == column ? static_cast< double >( row + 1 ) : 0.01 * static_cast< double >( 10 * ( row + 1 ) + column + 1 );
        }
    }
    rejectedSetBlockDataset.addObservationSetWithWeights( angular_position,
                                                          linkDefinition,
                                                          { ( Eigen::Vector2d( ) << 10.0, 11.0 ).finished( ),
                                                            ( Eigen::Vector2d( ) << 20.0, 21.0 ).finished( ),
                                                            ( Eigen::Vector2d( ) << 30.0, 31.0 ).finished( ) },
                                                          { 1.0, 2.0, 3.0 },
                                                          receiver,
                                                          ObservationWeightSettings::setBlock( fullSetWeightBlock ) );
    rejectedSetBlockDataset.rejectObservations( ObservationSelectionCondition< double, double >::timeBounds( 1.5, 2.5 ) );
    Eigen::MatrixXd expectedRejectedWeightBlock = Eigen::MatrixXd::Zero( 4, 4 );
    const std::vector< int > keptScalarIndices = { 0, 1, 4, 5 };
    for( unsigned int row = 0; row < keptScalarIndices.size( ); ++row )
    {
        for( unsigned int column = 0; column < keptScalarIndices.size( ); ++column )
        {
            expectedRejectedWeightBlock( row, column ) = fullSetWeightBlock( keptScalarIndices.at( row ), keptScalarIndices.at( column ) );
        }
    }

    // Projecting after rejection must keep only the active rows and columns of a full set weight block.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            rejectedSetBlockDataset.createEstimationFlattenedObservationData( ).getSparseWeightMatrix( ).toDense( ),
            expectedRejectedWeightBlock,
            1.0E-15 );
    rejectedSetBlockDataset.restoreObservations( ObservationSelectionCondition< double, double >::rejected( ) );

    // Restoring rejected rows must restore the full original weight block in the flattened data.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            rejectedSetBlockDataset.createEstimationFlattenedObservationData( ).getSparseWeightMatrix( ).toDense( ),
            fullSetWeightBlock,
            1.0E-15 );
}

/*!
 * Verifies the precedence between set-level, per-observation and extra block weights.
 *
 * Test outline: creates one vector-valued set with a full set-level block, then
 * overwrites part of it with an explicit per-observation block and finally
 * overwrites cross-observation entries with an extra scalar-component block. It
 * checks the resulting dense sparse matrix and confirms that conflicts report
 * the affected flattened/scalar-component indices.
 */
BOOST_AUTO_TEST_CASE( test_dataset_weight_precedence_and_conflict_warnings )
{
    const LinkDefinition linkDefinition = createOneWayLinkDefinition( "Station1" );

    ObservationDataset< double, double > dataset;
    const unsigned int setId = dataset.addObservationSet(
            angular_position,
            linkDefinition,
            { ( Eigen::Vector2d( ) << 10.0, 11.0 ).finished( ), ( Eigen::Vector2d( ) << 20.0, 21.0 ).finished( ) },
            { 1.0, 2.0 },
            receiver );
    const std::vector< unsigned int > observationIds = dataset.getObservationIdsForSet( setId );

    Eigen::MatrixXd setWeightBlock( 4, 4 );
    setWeightBlock << 1.0, 0.1, 0.2, 0.3, 0.1, 2.0, 0.4, 0.5, 0.2, 0.4, 3.0, 0.6, 0.3, 0.5, 0.6, 4.0;
    dataset.setWeightMatrixForSet( setId, setWeightBlock );

    Eigen::Matrix2d perObservationWeightBlock;
    perObservationWeightBlock << 10.0, 0.7, 0.8, 11.0;
    dataset.setWeightMatrixForObservation( observationIds.at( 0 ), perObservationWeightBlock );
    const Eigen::Vector2d perObservationDiagonalWeights = ( Eigen::Vector2d( ) << 12.0, 13.0 ).finished( );
    dataset.setConstantSingleObservationDiagonalWeight( ObservationSelectionCondition< double, double >::timeBounds( 1.5, 2.5 ),
                                                        perObservationDiagonalWeights );

    Eigen::Matrix2d extraWeightBlock;
    extraWeightBlock << 20.0, 21.0, 22.0, 23.0;
    dataset.setWeightBlock( { observationIds.at( 0 ) }, { observationIds.at( 1 ) }, extraWeightBlock );

    std::ostringstream warningStream;
    std::streambuf* originalWarningBuffer = std::cerr.rdbuf( warningStream.rdbuf( ) );
    const FlattenedObservationData< double, double > flattenedData = dataset.createEstimationFlattenedObservationData( );
    std::cerr.rdbuf( originalWarningBuffer );

    Eigen::MatrixXd expectedWeightMatrix = setWeightBlock;
    expectedWeightMatrix.block( 0, 0, 2, 2 ) = perObservationWeightBlock;
    expectedWeightMatrix.block( 2, 2, 2, 2 ) = perObservationDiagonalWeights.asDiagonal( );
    expectedWeightMatrix.block( 0, 2, 2, 2 ) = extraWeightBlock;
    expectedWeightMatrix.block( 2, 0, 2, 2 ) = extraWeightBlock.transpose( );

    // Set-level blocks provide the baseline, explicit matrix/vector observation weights overwrite that baseline,
    // and extra scalar-component blocks overwrite both lower-priority layers while preserving matrix symmetry.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( flattenedData.getSparseWeightMatrix( ).toDense( ), expectedWeightMatrix, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( flattenedData.getWeightVector( ), expectedWeightMatrix.diagonal( ), 1.0E-15 );

    const std::string warningText = warningStream.str( );

    // The per-observation layer conflicts with the set-level diagonal entry at the first scalar component.
    BOOST_CHECK( warningText.find( "flattened matrix row 0, column 0 (scalar component ids 0, 0)" ) != std::string::npos );

    // The per-observation layer also conflicts with an off-diagonal entry inside the first observation block.
    BOOST_CHECK( warningText.find( "flattened matrix row 0, column 1 (scalar component ids 0, 1)" ) != std::string::npos );

    // The vector observation-level update clears the within-observation off-diagonal set-level entry.
    BOOST_CHECK( warningText.find( "flattened matrix row 2, column 3 (scalar component ids 2, 3)" ) != std::string::npos );

    // The extra scalar-component block conflicts with the set-level block across the two observation rows.
    BOOST_CHECK( warningText.find( "flattened matrix row 0, column 2 (scalar component ids 0, 2)" ) != std::string::npos );

    // The warning must document the resolved precedence so users can understand which value remains active.
    BOOST_CHECK( warningText.find( "Precedence is set-level weights, then explicit per-observation weights, "
                                   "then extra scalar-component blocks" ) != std::string::npos );
}

/*!
 * Verifies that default weights introduced while appending rows remain implicit.
 *
 * Test outline: appends rows before setting a full set-level weight block. The
 * first dataset omits append weights, so the later set-level block should remain
 * authoritative for every row. The second dataset supplies append weights
 * explicitly, so only that appended row should override its diagonal entry.
 */
BOOST_AUTO_TEST_CASE( test_appended_default_weights_do_not_override_later_set_block )
{
    const LinkDefinition linkDefinition = createOneWayLinkDefinition( "Station1" );

    ObservationDataset< double, double > implicitAppendDataset;
    const unsigned int implicitSetId =
            implicitAppendDataset.addObservationSet( one_way_range,
                                                     linkDefinition,
                                                     { Eigen::Vector1d::Constant( 10.0 ), Eigen::Vector1d::Constant( 20.0 ) },
                                                     { 1.0, 2.0 },
                                                     receiver );
    implicitAppendDataset.addObservationsToSet( implicitSetId, { Eigen::Vector1d::Constant( 30.0 ) }, { 3.0 }, {}, {}, {}, true );

    Eigen::Matrix3d setWeightBlock;
    setWeightBlock << 1.0, 0.1, 0.2, 0.1, 2.0, 0.3, 0.2, 0.3, 3.0;
    implicitAppendDataset.setWeightMatrixForSet( implicitSetId, setWeightBlock );

    std::ostringstream implicitWarningStream;
    std::streambuf* originalWarningBuffer = std::cerr.rdbuf( implicitWarningStream.rdbuf( ) );
    const Eigen::MatrixXd implicitFlattenedWeightMatrix =
            implicitAppendDataset.createEstimationFlattenedObservationData( ).getSparseWeightMatrix( ).toDense( );
    std::cerr.rdbuf( originalWarningBuffer );

    // Appended rows without user-supplied weights should keep implicit unit defaults, so the set block should be unchanged.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( implicitFlattenedWeightMatrix, setWeightBlock, 1.0E-15 );
    // No conflict warning should be printed when implicit defaults are correctly skipped under a set-level block.
    BOOST_CHECK( implicitWarningStream.str( ).empty( ) );

    ObservationDataset< double, double > explicitAppendDataset;
    const unsigned int explicitSetId =
            explicitAppendDataset.addObservationSet( one_way_range,
                                                     linkDefinition,
                                                     { Eigen::Vector1d::Constant( 10.0 ), Eigen::Vector1d::Constant( 20.0 ) },
                                                     { 1.0, 2.0 },
                                                     receiver );
    explicitAppendDataset.addObservationsToSet(
            explicitSetId, { Eigen::Vector1d::Constant( 30.0 ) }, { 3.0 }, {}, { Eigen::Vector1d::Constant( 9.0 ) }, {}, true );
    explicitAppendDataset.setWeightMatrixForSet( explicitSetId, setWeightBlock );

    std::ostringstream explicitWarningStream;
    originalWarningBuffer = std::cerr.rdbuf( explicitWarningStream.rdbuf( ) );
    const Eigen::MatrixXd explicitFlattenedWeightMatrix =
            explicitAppendDataset.createEstimationFlattenedObservationData( ).getSparseWeightMatrix( ).toDense( );
    std::cerr.rdbuf( originalWarningBuffer );

    Eigen::Matrix3d expectedExplicitWeightMatrix = setWeightBlock;
    expectedExplicitWeightMatrix( 2, 2 ) = 9.0;

    // Explicit append weights should override the diagonal entry for the appended observation only.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( explicitFlattenedWeightMatrix, expectedExplicitWeightMatrix, 1.0E-15 );
    // The override should be reported with the flattened/scalar-component index that changed.
    BOOST_CHECK( explicitWarningStream.str( ).find( "flattened matrix row 2, column 2 (scalar component ids 2, 2)" ) != std::string::npos );
}

/*!
 * Verifies that row status and advanced weights survive dataset rebuilds.
 *
 * Test outline: creates rejected rows and off-diagonal weights, then exercises
 * keep/drop, explicit removal and append operations. It checks that rejecting a
 * row never deletes it, that active flattened data exclude it, and that compact,
 * per-observation, set-level and arbitrary weight blocks are preserved or
 * subsetted when the operation has a well-defined mapping.
 */
BOOST_AUTO_TEST_CASE( test_dataset_rebuild_preserves_status_and_weights )
{
    const LinkDefinition linkDefinition = createOneWayLinkDefinition( "Station1" );
    ObservationDataset< double, double > dataset;
    const unsigned int setId = dataset.addObservationSet(
            one_way_range,
            linkDefinition,
            { Eigen::Vector1d::Constant( 10.0 ), Eigen::Vector1d::Constant( 20.0 ), Eigen::Vector1d::Constant( 30.0 ) },
            { 1.0, 2.0, 3.0 },
            receiver,
            std::vector< Eigen::VectorXd >( ),
            nullptr,
            nullptr,
            { Eigen::Vector1d::Constant( 2.0 ), Eigen::Vector1d::Constant( 3.0 ), Eigen::Vector1d::Constant( 4.0 ) } );
    const std::vector< unsigned int > observationIds = dataset.getObservationIdsForSet( setId );
    dataset.setWeightMatrixForObservation( observationIds.at( 0 ), ( Eigen::MatrixXd( 1, 1 ) << 5.0 ).finished( ) );
    dataset.setWeightBlock( { observationIds.at( 0 ), observationIds.at( 2 ) },
                            { observationIds.at( 0 ), observationIds.at( 2 ) },
                            ( Eigen::Matrix2d( ) << 6.0, 0.25, 0.25, 7.0 ).finished( ),
                            {},
                            {} );
    dataset.rejectObservations( ObservationSelectionCondition< double, double >::timeBounds( 1.5, 2.5 ), "manual rejection" );

    // Rejected observations remain stored for computation/inspection, while the default estimation flattened data excludes them.
    BOOST_CHECK_EQUAL( dataset.getNumberOfObservations( ), 3 );
    BOOST_CHECK_EQUAL( dataset.createComputationFlattenedObservationData( true ).getObservationVector( ).size( ), 3 );
    BOOST_CHECK_EQUAL( dataset.createEstimationFlattenedObservationData( false ).getObservationVector( ).size( ), 2 );

    const std::shared_ptr< ObservationDataset< double, double > > keepAll =
            dataset.createNewAndKeep( ObservationSelectionCondition< double, double >::all( ) );

    // Copy-style reduction with all rows must preserve row rejection state and the complete sparse weight flattened data.
    BOOST_CHECK( !keepAll->getObservationRow( 1 ).isActive_ );
    BOOST_CHECK_EQUAL( keepAll->getObservationRow( 1 ).rejectionReason_, "manual rejection" );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( keepAll->createComputationFlattenedObservationData( true ).getSparseWeightMatrix( ).toDense( ),
                                       dataset.createComputationFlattenedObservationData( true ).getSparseWeightMatrix( ).toDense( ),
                                       1.0E-15 );

    const std::shared_ptr< ObservationDataset< double, double > > activeOnly =
            dataset.createNewAndDrop( ObservationSelectionCondition< double, double >::rejected( ) );

    // Explicit drop is the only operation here that physically removes the rejected observation from the new dataset.
    BOOST_CHECK_EQUAL( dataset.getNumberOfObservations( ), 3 );
    BOOST_CHECK_EQUAL( activeOnly->getNumberOfObservations( ), 2 );
    BOOST_CHECK_EQUAL( activeOnly->createEstimationFlattenedObservationData( true ).getObservationVector( ).size( ), 2 );

    dataset.removeObservationsFromSet( setId, std::vector< unsigned int >( { 0 } ) );

    // Removing a different row must not reset the rejection state of the remaining rejected observation.
    BOOST_CHECK_EQUAL( dataset.getNumberOfObservations( ), 2 );
    BOOST_CHECK( !dataset.getObservationRow( 0 ).isActive_ );
    BOOST_CHECK_EQUAL( dataset.getObservationRow( 0 ).rejectionReason_, "manual rejection" );
    BOOST_CHECK_EQUAL( dataset.createEstimationFlattenedObservationData( false ).getObservationVector( ).size( ), 1 );
    BOOST_CHECK_EQUAL( dataset.createComputationFlattenedObservationData( true ).getObservationVector( ).size( ), 2 );

    dataset.addObservationsToSet( setId,
                                  { Eigen::Vector1d::Constant( 40.0 ) },
                                  { 4.0 },
                                  std::vector< Eigen::VectorXd >( ),
                                  { Eigen::Vector1d::Constant( 8.0 ) },
                                  std::vector< Eigen::VectorXd >( ),
                                  true );

    // Appending a new row must preserve the old rejected row and assign the provided compact weight to the new row.
    BOOST_CHECK( !dataset.getObservationRow( 0 ).isActive_ );
    BOOST_CHECK_EQUAL( dataset.getObservationRow( 0 ).rejectionReason_, "manual rejection" );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            dataset.getWeightValue( dataset.getObservationIdsForSet( setId ).back( ) ), Eigen::Vector1d::Constant( 8.0 ), 1.0E-15 );

    ObservationDataset< double, double > setBlockDataset;
    const Eigen::Matrix3d setWeightBlock = ( Eigen::Matrix3d( ) << 1.0, 0.1, 0.2, 0.1, 2.0, 0.3, 0.2, 0.3, 3.0 ).finished( );
    const int setBlockSetId = setBlockDataset.addObservationSetWithWeights(
            one_way_range,
            linkDefinition,
            { Eigen::Vector1d::Constant( 10.0 ), Eigen::Vector1d::Constant( 20.0 ), Eigen::Vector1d::Constant( 30.0 ) },
            { 1.0, 2.0, 3.0 },
            receiver,
            ObservationWeightSettings::setBlock( setWeightBlock ) );

    const std::shared_ptr< ObservationDataset< double, double > > setBlockDatasetPointer =
            std::make_shared< ObservationDataset< double, double > >( setBlockDataset );
    const std::vector< std::shared_ptr< SingleObservationSet< double, double > > > splitSets =
            splitObservationSet( createSingleObservationSet( setBlockDatasetPointer, setBlockSetId ),
                                 observationSetSplitter( nb_observations_splitter, 2, 1 ),
                                 false );
    BOOST_REQUIRE_EQUAL( splitSets.size( ), 2 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            splitSets.at( 0 )->getObservationDataset( )->createEstimationFlattenedObservationData( ).getSparseWeightMatrix( ).toDense( ),
            setWeightBlock.topLeftCorner( 2, 2 ),
            1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            splitSets.at( 1 )->getObservationDataset( )->createEstimationFlattenedObservationData( ).getSparseWeightMatrix( ).toDense( ),
            Eigen::MatrixXd::Constant( 1, 1, setWeightBlock( 2, 2 ) ),
            1.0E-15 );

    const std::shared_ptr< ObservationCollection< double, double > > copiedCollection = createNewObservationCollection(
            createObservationCollection( std::make_shared< ObservationDataset< double, double > >( setBlockDataset ) ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            copiedCollection->getObservationDataset( )->createEstimationFlattenedObservationData( ).getSparseWeightMatrix( ).toDense( ),
            setWeightBlock,
            1.0E-15 );

    setBlockDataset.removeObservationsFromSet( setBlockSetId, std::vector< unsigned int >( { 1 } ) );
    const Eigen::Matrix2d expectedSubsetSetWeightBlock = ( Eigen::Matrix2d( ) << 1.0, 0.2, 0.2, 3.0 ).finished( );

    // Removing from a set-level block has a clear subsetting rule and must keep the corresponding rows/columns.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( setBlockDataset.createEstimationFlattenedObservationData( ).getSparseWeightMatrix( ).toDense( ),
                                       expectedSubsetSetWeightBlock,
                                       1.0E-15 );

    // Appending to a full set-level block has no unique correlation extension and must be rejected explicitly.
    BOOST_CHECK_THROW( setBlockDataset.addObservationsToSet( setBlockSetId, { Eigen::Vector1d::Constant( 40.0 ) }, { 4.0 } ),
                       std::runtime_error );

    ObservationDataset< double, double > moveSourceDataset;
    const unsigned int moveSourceSetId = moveSourceDataset.addObservationSet(
            one_way_range,
            linkDefinition,
            { Eigen::Vector1d::Constant( 10.0 ), Eigen::Vector1d::Constant( 20.0 ), Eigen::Vector1d::Constant( 30.0 ) },
            { 1.0, 2.0, 3.0 },
            receiver );
    const std::vector< unsigned int > moveSourceObservationIds = moveSourceDataset.getObservationIdsForSet( moveSourceSetId );
    const Eigen::Matrix2d movedWeightBlock = ( Eigen::Matrix2d( ) << 5.0, 0.4, 0.4, 7.0 ).finished( );
    moveSourceDataset.setWeightBlock( { moveSourceObservationIds.at( 0 ), moveSourceObservationIds.at( 2 ) },
                                      { moveSourceObservationIds.at( 0 ), moveSourceObservationIds.at( 2 ) },
                                      movedWeightBlock );
    ObservationDataset< double, double > moveTargetDataset;
    const unsigned int moveTargetSetId = moveTargetDataset.addObservationSet(
            one_way_range, linkDefinition, std::vector< Eigen::VectorXd >( ), std::vector< double >( ), receiver );
    moveSourceDataset.moveObservationsToSet( moveSourceSetId, moveTargetDataset, moveTargetSetId, { 0, 2 }, false );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            moveTargetDataset.createEstimationFlattenedObservationData( ).getSparseWeightMatrix( ).toDense( ), movedWeightBlock, 1.0E-15 );
}

/*!
 * Verifies that copying a set between datasets preserves row status and sparse extra weights.
 *
 * Test outline: creates a source set with a rejected observation and an
 * off-diagonal extra weight block, copies it into a non-empty target dataset,
 * and verifies that the copied rows keep their active/rejected state and that
 * the extra weight block is remapped from source scalar-component ids to target
 * scalar-component ids.
 */
BOOST_AUTO_TEST_CASE( test_add_observation_set_from_dataset_preserves_status_and_extra_weights )
{
    const LinkDefinition linkDefinition = createOneWayLinkDefinition( "Station1" );

    ObservationDataset< double, double > sourceDataset;
    const unsigned int sourceSetId = sourceDataset.addObservationSet( angular_position,
                                                                      linkDefinition,
                                                                      { ( Eigen::Vector2d( ) << 10.0, 11.0 ).finished( ),
                                                                        ( Eigen::Vector2d( ) << 20.0, 21.0 ).finished( ),
                                                                        ( Eigen::Vector2d( ) << 30.0, 31.0 ).finished( ) },
                                                                      { 1.0, 2.0, 3.0 },
                                                                      receiver );
    const std::vector< unsigned int > sourceObservationIds = sourceDataset.getObservationIdsForSet( sourceSetId );

    const Eigen::Matrix2d sourceExtraWeightBlock = ( Eigen::Matrix2d( ) << 4.0, 0.2, 0.3, 5.0 ).finished( );
    sourceDataset.setWeightBlock( { sourceObservationIds.at( 0 ) }, { sourceObservationIds.at( 2 ) }, sourceExtraWeightBlock );
    sourceDataset.rejectObservations( ObservationSelectionCondition< double, double >::timeBounds( 1.5, 2.5 ), "copy rejection" );

    // The source setup must contain exactly the state that addObservationSetFromDataset is expected to preserve.
    BOOST_CHECK( !sourceDataset.getObservationRow( sourceObservationIds.at( 1 ) ).isActive_ );
    BOOST_CHECK_EQUAL( sourceDataset.getObservationRow( sourceObservationIds.at( 1 ) ).rejectionReason_, "copy rejection" );
    BOOST_REQUIRE_EQUAL( sourceDataset.getExtraWeightBlocks( ).size( ), 2 );

    ObservationDataset< double, double > targetDataset;
    targetDataset.addObservationSet( one_way_range, linkDefinition, { Eigen::Vector1d::Constant( 1.0 ) }, { 0.0 }, receiver );
    const unsigned int copiedSetId = targetDataset.addObservationSetFromDataset( sourceDataset, sourceSetId );
    const std::vector< unsigned int > copiedObservationIds = targetDataset.getObservationIdsForSet( copiedSetId );

    // The target starts non-empty so copied observation and scalar-component ids differ from the source ids.
    BOOST_CHECK_EQUAL( copiedSetId, 1 );
    BOOST_CHECK( copiedObservationIds.at( 0 ) != sourceObservationIds.at( 0 ) );
    BOOST_CHECK( targetDataset.getObservationRow( copiedObservationIds.at( 0 ) ).firstScalarComponent_ !=
                 sourceDataset.getObservationRow( sourceObservationIds.at( 0 ) ).firstScalarComponent_ );

    // Active rows must remain active after copying, not just retain their observation values and times.
    BOOST_CHECK( targetDataset.getObservationRow( copiedObservationIds.at( 0 ) ).isActive_ );
    BOOST_CHECK( targetDataset.getObservationRow( copiedObservationIds.at( 2 ) ).isActive_ );

    // Rejected rows must remain stored but inactive, with the original rejection reason preserved for diagnostics.
    BOOST_CHECK( !targetDataset.getObservationRow( copiedObservationIds.at( 1 ) ).isActive_ );
    BOOST_CHECK_EQUAL( targetDataset.getObservationRow( copiedObservationIds.at( 1 ) ).rejectionReason_, "copy rejection" );
    checkIds( targetDataset.getObservationIdsMatchingCondition( ObservationSelectionCondition< double, double >::rejected( ) ),
              { copiedObservationIds.at( 1 ) } );

    // Default estimation data must still exclude the copied rejected row while keeping the unrelated pre-existing target row.
    BOOST_CHECK_EQUAL( targetDataset.getNumberOfObservations( ), 4 );
    BOOST_CHECK_EQUAL( targetDataset.createEstimationFlattenedObservationData( false ).getObservationVector( ).size( ), 5 );
    BOOST_CHECK_EQUAL( targetDataset.createComputationFlattenedObservationData( true ).getObservationVector( ).size( ), 7 );

    BOOST_REQUIRE_EQUAL( targetDataset.getExtraWeightBlocks( ).size( ), 2 );
    const ObservationWeightBlock& copiedWeightBlock = targetDataset.getExtraWeightBlocks( ).at( 0 );
    const ObservationWeightBlock& copiedTransposedWeightBlock = targetDataset.getExtraWeightBlocks( ).at( 1 );
    const ObservationDatasetRow< double >& copiedFirstRow = targetDataset.getObservationRow( copiedObservationIds.at( 0 ) );
    const ObservationDatasetRow< double >& copiedThirdRow = targetDataset.getObservationRow( copiedObservationIds.at( 2 ) );

    // The stored sparse blocks must be remapped to the copied rows' scalar-component ids, not left pointing at source ids.
    checkIds( copiedWeightBlock.rowScalarComponentIds_,
              { copiedFirstRow.firstScalarComponent_, copiedFirstRow.firstScalarComponent_ + 1 } );
    checkIds( copiedWeightBlock.columnScalarComponentIds_,
              { copiedThirdRow.firstScalarComponent_, copiedThirdRow.firstScalarComponent_ + 1 } );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( copiedWeightBlock.weightBlock_, sourceExtraWeightBlock, 1.0E-15 );
    checkIds( copiedTransposedWeightBlock.rowScalarComponentIds_,
              { copiedThirdRow.firstScalarComponent_, copiedThirdRow.firstScalarComponent_ + 1 } );
    checkIds( copiedTransposedWeightBlock.columnScalarComponentIds_,
              { copiedFirstRow.firstScalarComponent_, copiedFirstRow.firstScalarComponent_ + 1 } );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( copiedTransposedWeightBlock.weightBlock_, sourceExtraWeightBlock.transpose( ), 1.0E-15 );

    const FlattenedObservationData< double, double > copiedFlattenedData = targetDataset.createComputationFlattenedObservationData( true );
    const Eigen::MatrixXd copiedDenseWeightMatrix = copiedFlattenedData.getSparseWeightMatrix( ).toDense( );
    const int copiedFirstRowStart = copiedFlattenedData.getFlattenedRow( copiedObservationIds.at( 0 ), 0 );
    const int copiedThirdRowStart = copiedFlattenedData.getFlattenedRow( copiedObservationIds.at( 2 ), 0 );

    // Downstream consumers use the flattened sparse matrix, so the remapped block must materialize at the copied row positions.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            copiedDenseWeightMatrix.block( copiedFirstRowStart, copiedThirdRowStart, 2, 2 ), sourceExtraWeightBlock, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            copiedDenseWeightMatrix.block( copiedThirdRowStart, copiedFirstRowStart, 2, 2 ), sourceExtraWeightBlock.transpose( ), 1.0E-15 );
}

BOOST_AUTO_TEST_CASE( test_dataset_copies_own_mutable_metadata_and_support_self_copy )
{
    const LinkDefinition linkDefinition = createOneWayLinkDefinition( "Station1" );
    const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > elevationSettings =
            simulation_setup::elevationAngleDependentVariable( receiver, LinkEndId( "Vehicle", "" ) );
    const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping > bookkeeping =
            std::make_shared< simulation_setup::ObservationDependentVariableBookkeeping >( one_way_range, linkDefinition );
    bookkeeping->addDependentVariable( elevationSettings );
    const std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings =
            std::make_shared< ObservationAncillarySimulationSettings >( );
    ancillarySettings->setAncillaryDoubleData( doppler_integration_time, 10.0 );

    ObservationDataset< double, double > sourceDataset;
    const unsigned int sourceSetId = sourceDataset.addObservationSet(
            one_way_range, linkDefinition, { Eigen::Vector1d::Constant( 1.0 ) }, { 1.0 }, receiver, {}, bookkeeping, ancillarySettings );
    const unsigned int angularSetId =
            sourceDataset.addObservationSet( angular_position, linkDefinition, { Eigen::Vector2d( 2.0, 3.0 ) }, { 2.0 }, receiver );

    ObservationDataset< double, double > copiedDataset;
    const unsigned int copiedSetId = copiedDataset.addObservationSetFromDataset( sourceDataset, sourceSetId );
    copiedDataset.getDependentVariableBookkeeping( copiedDataset.getObservationSetMetadata( copiedSetId ).dependentVariableLayoutId_ )
            ->clearSettings( );
    copiedDataset.getAncillarySettings( copiedDataset.getObservationSetMetadata( copiedSetId ).ancillarySettingsId_ )
            ->setAncillaryDoubleData( doppler_integration_time, 20.0 );

    BOOST_CHECK_EQUAL(
            sourceDataset
                    .getDependentVariableBookkeeping( sourceDataset.getObservationSetMetadata( sourceSetId ).dependentVariableLayoutId_ )
                    ->getTotalDependentVariableSize( ),
            1 );
    BOOST_CHECK_EQUAL( sourceDataset.getAncillarySettings( sourceDataset.getObservationSetMetadata( sourceSetId ).ancillarySettingsId_ )
                               ->getAncillaryDoubleData( doppler_integration_time ),
                       10.0 );

    const std::shared_ptr< ObservationDataset< double, double > > reducedDataset =
            sourceDataset.createNewAndKeep( ObservationSelectionCondition< double, double >::observableType( one_way_range ) );
    reducedDataset->getDependentVariableBookkeeping( reducedDataset->getObservationSetMetadata( 0 ).dependentVariableLayoutId_ )
            ->clearSettings( );
    BOOST_CHECK_EQUAL(
            sourceDataset
                    .getDependentVariableBookkeeping( sourceDataset.getObservationSetMetadata( sourceSetId ).dependentVariableLayoutId_ )
                    ->getTotalDependentVariableSize( ),
            1 );

    const std::shared_ptr< ObservationDataset< double, double > > singleSetDataset = createObservationDataset(
            createSingleObservationSet( std::make_shared< ObservationDataset< double, double > >( sourceDataset ), angularSetId ) );
    BOOST_CHECK_EQUAL( singleSetDataset->getNumberOfObservationSets( ), 1 );
    BOOST_CHECK_EQUAL( singleSetDataset->getObservationSetMetadata( 0 ).observableType_, angular_position );

    ObservationDataset< double, double > copyOnWriteDataset;
    copyOnWriteDataset.addObservationSet(
            one_way_range, linkDefinition, { Eigen::Vector1d::Constant( 1.0 ) }, { 1.0 }, receiver, {}, bookkeeping );
    copyOnWriteDataset.addObservationSet(
            one_way_range, linkDefinition, { Eigen::Vector1d::Constant( 2.0 ) }, { 2.0 }, receiver, {}, bookkeeping );
    const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > transmitterElevationSettings =
            simulation_setup::elevationAngleDependentVariable( transmitter, LinkEndId( "Earth", "Station1" ) );
    copyOnWriteDataset.addDependentVariableToSets( transmitterElevationSettings,
                                                   ObservationSelectionCondition< double, double >::timeBounds( 0.5, 1.5 ) );

    BOOST_CHECK_EQUAL(
            copyOnWriteDataset
                    .getDependentVariableBookkeeping( copyOnWriteDataset.getObservationSetMetadata( 0 ).dependentVariableLayoutId_ )
                    ->getTotalDependentVariableSize( ),
            2 );
    BOOST_CHECK_EQUAL(
            copyOnWriteDataset
                    .getDependentVariableBookkeeping( copyOnWriteDataset.getObservationSetMetadata( 1 ).dependentVariableLayoutId_ )
                    ->getTotalDependentVariableSize( ),
            1 );
    BOOST_CHECK_EQUAL( bookkeeping->getTotalDependentVariableSize( ), 1 );

    ObservationDataset< double, double > selfCopyDataset;
    const unsigned int selfCopySetId = selfCopyDataset.addObservationSet(
            one_way_range, linkDefinition, { Eigen::Vector1d::Constant( 1.0 ), Eigen::Vector1d::Constant( 2.0 ) }, { 1.0, 2.0 }, receiver );
    const std::vector< unsigned int > selfCopyObservationIds = selfCopyDataset.getObservationIdsForSet( selfCopySetId );
    selfCopyDataset.setWeightBlock(
            selfCopyObservationIds, selfCopyObservationIds, ( Eigen::Matrix2d( ) << 3.0, 0.25, 0.25, 4.0 ).finished( ) );
    const unsigned int newSelfCopiedSetId = selfCopyDataset.addObservationSetFromDataset( selfCopyDataset, selfCopySetId );

    BOOST_CHECK_EQUAL( newSelfCopiedSetId, 1 );
    BOOST_CHECK_EQUAL( selfCopyDataset.getNumberOfObservationSets( ), 2 );
    BOOST_CHECK_EQUAL( selfCopyDataset.getExtraWeightBlocks( ).size( ), 2 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( selfCopyDataset.getObservationVectorForSet( newSelfCopiedSetId ),
                                       selfCopyDataset.getObservationVectorForSet( selfCopySetId ),
                                       1.0E-15 );
}

/*!
 * Verifies the sparse-weight least-squares overload.
 *
 * Test outline: solves a small weighted least-squares problem with off-diagonal
 * sparse weights and compares the returned parameter update and normal matrix
 * against an independently assembled dense reference calculation.
 */
BOOST_AUTO_TEST_CASE( test_sparse_weighted_least_squares )
{
    Eigen::MatrixXd designMatrix( 3, 2 );
    designMatrix << 1.0, 2.0, 3.0, -1.0, 0.5, 4.0;
    Eigen::VectorXd residuals( 3 );
    residuals << 0.2, -0.4, 0.7;

    Eigen::SparseMatrix< double > sparseWeights( 3, 3 );
    std::vector< Eigen::Triplet< double > > weightTriplets;
    weightTriplets.emplace_back( 0, 0, 2.0 );
    weightTriplets.emplace_back( 1, 1, 3.0 );
    weightTriplets.emplace_back( 2, 2, 4.0 );
    weightTriplets.emplace_back( 0, 2, 0.3 );
    weightTriplets.emplace_back( 2, 0, 0.3 );
    sparseWeights.setFromTriplets( weightTriplets.begin( ), weightTriplets.end( ) );

    const Eigen::MatrixXd denseWeights = sparseWeights.toDense( );
    const Eigen::MatrixXd expectedNormalMatrix = designMatrix.transpose( ) * denseWeights * designMatrix;
    const Eigen::VectorXd expectedRightHandSide = designMatrix.transpose( ) * denseWeights * residuals;
    const Eigen::VectorXd expectedParameterUpdate =
            linear_algebra::solveSystemOfEquationsWithSvd( expectedNormalMatrix, expectedRightHandSide );

    const std::pair< Eigen::VectorXd, Eigen::MatrixXd > sparseLeastSquaresOutput =
            linear_algebra::performLeastSquaresAdjustmentFromDesignMatrix( designMatrix, residuals, sparseWeights );

    // Sparse-weight least squares must produce the same normal matrix and parameter update as the dense reference.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( sparseLeastSquaresOutput.second, expectedNormalMatrix, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( sparseLeastSquaresOutput.first, expectedParameterUpdate, 1.0E-15 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
