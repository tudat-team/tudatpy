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

#include <algorithm>
#include <cmath>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/orbit_determination/podInputOutputTypes.h"
#include "tudat/math/basic/leastSquaresEstimation.h"
#include "tudat/simulation/estimation_setup/observationOutputSettings.h"
#include "tudat/simulation/estimation_setup/observationDataset.h"
#include "tudat/simulation/estimation_setup/podProcessing.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::observation_models;

//! Test observation-dataset storage, selection, mutation, weighting, and inspection.
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

Eigen::VectorXd concatenateInspectionValues( const std::vector< Eigen::VectorXd >& values )
{
    Eigen::Index size = 0;
    for( const auto& value : values )
    {
        size += value.size( );
    }
    Eigen::VectorXd result( size );
    Eigen::Index start = 0;
    for( const auto& value : values )
    {
        result.segment( start, value.size( ) ) = value;
        start += value.size( );
    }
    return result;
}

/*!
 * Verifies the primary ObservationDataset storage model and observation vector data output.
 *
 * Test outline: creates range, angular-position and position sets with explicit
 * weights and residuals. It checks set ids, row/scalar-component bookkeeping,
 * observation vector ordering, time replication for vector
 * observables and weight vector assembly.
 */
//! Verify dataset storage, scalar-vector projection, and weights.
BOOST_AUTO_TEST_CASE( test_dataset_storage_and_observation_vector_data )
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
    BOOST_CHECK_EQUAL( dataset.getObservationRow( 3 ).firstScalarComponent_, 4 );
    BOOST_CHECK_EQUAL( dataset.getObservationRow( 3 ).scalarSize_, 2 );
    BOOST_CHECK_EQUAL( dataset.getScalarComponentRow( 5 ).observationId_, 3 );
    BOOST_CHECK_EQUAL( dataset.getScalarComponentRow( 5 ).componentIndex_, 1 );

    Eigen::VectorXd expectedObservations( 9 );
    expectedObservations << 10.0, 11.0, 20.0, 21.0, 22.0, 23.0, 30.0, 31.0, 32.0;
    const ObservationVectorData< double, double > observationVectorData = dataset.createObservationVectorData( );

    // These sets were inserted in the established observable/link order.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( observationVectorData.getObservationVector( ), expectedObservations, 1.0E-15 );

    const std::vector< double > expectedTimes = { 1.0, 2.0, 3.0, 3.0, 4.0, 4.0, 5.0, 5.0, 5.0 };
    const std::vector< double >& observationVectorDataTimes = observationVectorData.getTimes( );

    // Vector-valued observables must repeat their event time once per scalar component.
    BOOST_CHECK_EQUAL_COLLECTIONS(
            observationVectorDataTimes.begin( ), observationVectorDataTimes.end( ), expectedTimes.begin( ), expectedTimes.end( ) );

    // Vector-data component positions come from the same mapping as the vectors.
    const std::vector< int > expectedStarts = { 0, 1, 2, 4, 6 };
    const std::vector< int > expectedSizes = { 1, 1, 2, 2, 3 };
    for( unsigned int id = 0; id < expectedStarts.size( ); ++id )
    {
        for( int component = 0; component < expectedSizes.at( id ); ++component )
        {
            BOOST_CHECK_EQUAL( observationVectorData.getVectorRow( id, component ), expectedStarts.at( id ) + component );
        }
    }

    Eigen::VectorXd expectedWeights( 9 );
    expectedWeights << 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0;

    // Compact weights must be assembled in the same scalar-component order as the observations.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( observationVectorData.getWeightVector( ), expectedWeights, 1.0E-15 );
}


/*!
 * Verifies empty-set behavior.
 *
 * Test outline: creates an empty observation set and checks that all counters
 * and observation vector data are well-defined, that the set survives simulation
 * and residual computation, and that row-condition link updates leave it unchanged.
 */
//! Verify empty observation sets remain valid metadata containers.
BOOST_AUTO_TEST_CASE( test_dataset_empty_sets )
{
    const LinkDefinition linkDefinition = createOneWayLinkDefinition( "Station1" );
    ObservationDataset< double, double > dataset;

    const int emptySetId = dataset.addObservationSet(
            one_way_range, linkDefinition, std::vector< Eigen::VectorXd >( ), std::vector< double >( ), receiver );

    // Empty sets must remain valid metadata containers with zero-length observation vector data.
    BOOST_CHECK_EQUAL( dataset.getNumberOfObservationSets( ), 1 );
    BOOST_CHECK_EQUAL( dataset.getNumberOfObservationsForSet( emptySetId ), 0 );
    BOOST_CHECK_EQUAL( dataset.getTotalScalarSizeForSet( emptySetId ), 0 );
    BOOST_CHECK_EQUAL( dataset.createObservationVectorData( ).getObservationVector( ).size( ), 0 );

    // An empty metadata group needs neither a simulator nor an environment,
    // and must retain its set position through simulation and residual update.
    auto emptyDataset = std::make_shared< ObservationDataset<> >( dataset );
    const simulation_setup::SystemOfBodies emptyBodies;
    const auto emptySettings = simulation_setup::getObservationSimulationSettingsFromObservationDataset( emptyDataset, emptyBodies );
    const auto simulatedEmpty = simulation_setup::simulateObservationDataset< double, double >( emptySettings, {}, emptyBodies );
    BOOST_CHECK_EQUAL( simulatedEmpty->getNumberOfObservationSets( ), 1 );
    BOOST_CHECK_EQUAL( simulatedEmpty->getNumberOfObservationsForSet( emptySetId ), 0 );
    BOOST_CHECK( simulatedEmpty->getLinkDefinition( simulatedEmpty->getObservationSetMetadata( emptySetId ).linkDefinitionId_ ) ==
                 linkDefinition );
    BOOST_CHECK_NO_THROW( simulation_setup::computeResidualsAndDependentVariables( emptyDataset, {}, emptyBodies ) );

    dataset.setLinkEndReferencePoint( "Earth", "RenamedStation", transmitter, ObservationSelectionCondition< double, double >::all( ) );
    BOOST_CHECK_EQUAL( dataset.getLinkDefinition( dataset.getObservationSetMetadata( emptySetId ).linkDefinitionId_ )
                               .linkEnds_.at( transmitter )
                               .getReferencePointName( ),
                       "Station1" );
}

/*!
 * Verifies sorting and duplicate removal.
 *
 * Test outline: adds intentionally unsorted and duplicate observations, removes
 * duplicates, and validates sorted row order and value/residual selection.
 */
//! Verify duplicate removal and row selection.
BOOST_AUTO_TEST_CASE( test_dataset_duplicate_removal_and_selection )
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
}

/*!
 * Verifies row-level ObservationSelectionCondition selection logic.
 *
 * Test outline: creates observations with link metadata, values, residuals, weights,
 * dependent variables and active/rejected status. It checks that each condition
 * selects exactly the intended observation ids and that combined conditions can
 * be used to build consistent snapshots.
 */

/*!
 * Verifies that snapshot ordered vector data follows dataset ordered-output order.
 *
 * Test outline: inserts angular-position rows before one-way range rows, so
 * dataset row order differs from ordered observation-vector-data order. Snapshot
 * selections retain dataset order while numerical vector data use the established
 * observable/link/set ordering, including after interleaved appends.
 */

/*!
 * Verifies snapshots, rejection/restoration and reduced dataset creation.
 *
 * Test outline: builds a dataset, selects rows with conditions, creates snapshots and
 * reduced datasets, rejects/restores observations and checks observation vector data sizes.
 * It also confirms that snapshots survive structural mutations.
 */


/*!
 * Verifies compact, block and sparse off-diagonal weight handling.
 *
 * Test outline: covers scalar weights, per-observation blocks, set-level blocks,
 * observation-id selected cross blocks, symmetric transpose insertion,
 * component selection and rejection/restoration of weighted observations. It
 * checks both compact weight vectors and materialized sparse observation vector data
 * matrices.
 */

/*!
 * Verifies the precedence between set-level, per-observation and extra block weights.
 *
 * Test outline: creates one vector-valued set with a full set-level block, then
 * overwrites part of it with an explicit per-observation block and finally
 * overwrites cross-observation entries with an extra scalar-component block. It
 * checks the resulting dense sparse matrix and confirms that conflicts report
 * the affected observation-vector and scalar-component indices.
 */

/*!
 * Verifies that later assignments consistently replace appended weights.
 *
 * Test outline: appends rows before setting a full set-level weight block. The
 * first dataset omits append weights, and the second supplies them explicitly.
 * A later set-level assignment must replace the addressed matrix in both cases.
 */

/*!
 * Verifies that row status and advanced weights survive dataset rebuilds.
 *
 * Test outline: creates rejected rows and off-diagonal weights, then exercises
 * keep/drop, explicit removal and append operations. It checks that rejecting a
 * row never deletes it, that active observation vector data exclude it, and that compact,
 * per-observation, set-level and arbitrary weight blocks are preserved or
 * subsetted when the operation has a well-defined mapping.
 */
//! Verify dataset rebuilds preserve row status and the restricted weight matrix.
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

    // Rejected observations remain stored for computation/inspection, while the default estimation observation vector data excludes them.
    BOOST_CHECK_EQUAL( dataset.getNumberOfObservations( ), 3 );
    BOOST_CHECK_EQUAL( dataset.createObservationVectorData( true ).getObservationVector( ).size( ), 3 );
    BOOST_CHECK_EQUAL( dataset.createObservationVectorData( false ).getObservationVector( ).size( ), 2 );

    const std::shared_ptr< ObservationDataset< double, double > > keepAll =
            dataset.createNewAndKeep( ObservationSelectionCondition< double, double >::all( ) );

    // Copy-style reduction with all rows must preserve row rejection state and the complete sparse weight observation vector data.
    BOOST_CHECK( !keepAll->getObservationRow( 1 ).isActive_ );
    BOOST_CHECK_EQUAL( keepAll->getObservationRow( 1 ).rejectionReason_, "manual rejection" );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( keepAll->createObservationVectorData( true ).getSparseWeightMatrix( ).toDense( ),
                                       dataset.createObservationVectorData( true ).getSparseWeightMatrix( ).toDense( ),
                                       1.0E-15 );

    const std::shared_ptr< ObservationDataset< double, double > > activeOnly =
            dataset.createNewAndDrop( ObservationSelectionCondition< double, double >::rejected( ) );

    // Explicit drop is the only operation here that physically removes the rejected observation from the new dataset.
    BOOST_CHECK_EQUAL( dataset.getNumberOfObservations( ), 3 );
    BOOST_CHECK_EQUAL( activeOnly->getNumberOfObservations( ), 2 );
    BOOST_CHECK_EQUAL( activeOnly->createObservationVectorData( true ).getObservationVector( ).size( ), 2 );

    dataset.removeObservations( ObservationSelectionCondition< double, double >::timeLessThan( 1.5 ) );

    // Removing a different row must not reset the rejection state of the remaining rejected observation.
    BOOST_CHECK_EQUAL( dataset.getNumberOfObservations( ), 2 );
    BOOST_CHECK( !dataset.getObservationRow( observationIds.at( 1 ) ).isActive_ );
    BOOST_CHECK_EQUAL( dataset.getObservationRow( observationIds.at( 1 ) ).rejectionReason_, "manual rejection" );
    BOOST_CHECK_EQUAL( dataset.createObservationVectorData( false ).getObservationVector( ).size( ), 1 );
    BOOST_CHECK_EQUAL( dataset.createObservationVectorData( true ).getObservationVector( ).size( ), 2 );

    dataset.addObservationsToSet( setId,
                                  { Eigen::Vector1d::Constant( 40.0 ) },
                                  { 4.0 },
                                  std::vector< Eigen::VectorXd >( ),
                                  { Eigen::Vector1d::Constant( 8.0 ) },
                                  std::vector< Eigen::VectorXd >( ),
                                  true );

    // Appending a new row must preserve the old rejected row and assign the provided compact weight to the new row.
    BOOST_CHECK( !dataset.getObservationRow( observationIds.at( 1 ) ).isActive_ );
    BOOST_CHECK_EQUAL( dataset.getObservationRow( observationIds.at( 1 ) ).rejectionReason_, "manual rejection" );
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
            setBlockWeights( setWeightBlock ) );

    setBlockDataset.removeObservations( ObservationSelectionCondition< double, double >::timeBounds( 1.5, 2.5 ) );
    const Eigen::Matrix2d expectedSubsetSetWeightBlock = ( Eigen::Matrix2d( ) << 1.0, 0.2, 0.2, 3.0 ).finished( );

    // Removing from a set-level block has a clear subsetting rule and must keep the corresponding rows/columns.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            setBlockDataset.createObservationVectorData( ).getSparseWeightMatrix( ).toDense( ), expectedSubsetSetWeightBlock, 1.0E-15 );

    // Appended rows start with unit diagonal weights and no correlations with existing rows.
    const auto survivorIds = setBlockDataset.getObservationIdsForSet( setBlockSetId );
    setBlockDataset.addObservationsToSet( setBlockSetId, { Eigen::Vector1d::Constant( 40.0 ) }, { 4.0 } );
    Eigen::Matrix3d expectedAppended = Eigen::Matrix3d::Identity( );
    expectedAppended.topLeftCorner< 2, 2 >( ) = expectedSubsetSetWeightBlock;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            setBlockDataset.createObservationVectorData( ).getSparseWeightMatrix( ).toDense( ), expectedAppended, 1.0E-15 );
    BOOST_CHECK_EQUAL( setBlockDataset.getObservationIdsForSet( setBlockSetId ).at( 0 ), survivorIds.at( 0 ) );
    BOOST_CHECK_EQUAL( setBlockDataset.getObservationIdsForSet( setBlockSetId ).at( 1 ), survivorIds.at( 1 ) );
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
//! Verify copying a set preserves rejection status and off-diagonal weights.
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
    BOOST_CHECK( sourceDataset.hasExtraWeightBlocks( ) );

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
    BOOST_CHECK_EQUAL( targetDataset.createObservationVectorData( false ).getObservationVector( ).size( ), 5 );
    BOOST_CHECK_EQUAL( targetDataset.createObservationVectorData( true ).getObservationVector( ).size( ), 7 );

    const ObservationVectorData< double, double > copiedObservationVectorData = targetDataset.createObservationVectorData( true );
    const Eigen::MatrixXd copiedDenseWeightMatrix = copiedObservationVectorData.getSparseWeightMatrix( ).toDense( );
    const int copiedFirstRowStart = copiedObservationVectorData.getVectorRow( copiedObservationIds.at( 0 ), 0 );
    const int copiedThirdRowStart = copiedObservationVectorData.getVectorRow( copiedObservationIds.at( 2 ), 0 );

    // Downstream consumers use the observation-vector sparse matrix, so the remapped block must materialize at the copied row positions.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            copiedDenseWeightMatrix.block( copiedFirstRowStart, copiedThirdRowStart, 2, 2 ), sourceExtraWeightBlock, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            copiedDenseWeightMatrix.block( copiedThirdRowStart, copiedFirstRowStart, 2, 2 ), sourceExtraWeightBlock.transpose( ), 1.0E-15 );
}

//! Verify dataset copies own mutable metadata and safely support self-copy operations.
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
    ObservationDataset< double, double > copyConstructedDataset( sourceDataset );
    copyConstructedDataset.getAncillarySettingsForSet( sourceSetId )->setAncillaryDoubleData( doppler_integration_time, 30.0 );

    BOOST_CHECK_EQUAL(
            sourceDataset
                    .getDependentVariableBookkeeping( sourceDataset.getObservationSetMetadata( sourceSetId ).dependentVariableLayoutId_ )
                    ->getTotalDependentVariableSize( ),
            1 );
    BOOST_CHECK_EQUAL( sourceDataset.getAncillarySettings( sourceDataset.getObservationSetMetadata( sourceSetId ).ancillarySettingsId_ )
                               ->getAncillaryDoubleData( doppler_integration_time ),
                       10.0 );

    // Filtered datasets must own copies of the bookkeeping and of the dependent-variable settings it contains.
    const std::shared_ptr< ObservationDataset< double, double > > reducedDataset =
            sourceDataset.createNewAndKeep( ObservationSelectionCondition< double, double >::observableType( one_way_range ) );
    const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& reducedBookkeeping =
            reducedDataset->getDependentVariableBookkeeping( reducedDataset->getObservationSetMetadata( 0 ).dependentVariableLayoutId_ );
    reducedBookkeeping->getDependentVariableSettings( ).front( )->linkEndId_ = LinkEndId( "Moon", "" );
    reducedBookkeeping->clearSettings( );
    const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& sourceBookkeeping =
            sourceDataset.getDependentVariableBookkeeping(
                    sourceDataset.getObservationSetMetadata( sourceSetId ).dependentVariableLayoutId_ );
    BOOST_CHECK_EQUAL( sourceBookkeeping->getTotalDependentVariableSize( ), 1 );
    BOOST_CHECK_EQUAL( sourceBookkeeping->getDependentVariableSettings( ).front( )->linkEndId_.bodyName_, "Vehicle" );

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
    Eigen::Matrix4d expectedCopiedWeights = Eigen::Matrix4d::Zero( );
    expectedCopiedWeights.topLeftCorner< 2, 2 >( ) << 3.0, 0.25, 0.25, 4.0;
    expectedCopiedWeights.bottomRightCorner< 2, 2 >( ) = expectedCopiedWeights.topLeftCorner< 2, 2 >( );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            selfCopyDataset.createObservationVectorData( ).getSparseWeightMatrix( ).toDense( ), expectedCopiedWeights, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( selfCopyDataset.getObservationVectorForSet( newSelfCopiedSetId ),
                                       selfCopyDataset.getObservationVectorForSet( selfCopySetId ),
                                       1.0E-15 );
}

/*!
 * Verifies the sparse-weight least-squares overload.
 *
 * Test outline: solves a small weighted least-squares problem with off-diagonal
 * sparse weights, with a dense and with a sparse design matrix, and compares the returned parameter update and normal
 * matrix, and the covariance contribution of a consider parameter, against an independently assembled dense reference.
 */

//! Verify removing rows does not renumber surviving observation identities.
BOOST_AUTO_TEST_CASE( test_removal_preserves_observation_identity )
{
    ObservationDataset<> dataset;
    const unsigned int setId = dataset.addObservationSet(
            one_way_range,
            createOneWayLinkDefinition( "Station1" ),
            { Eigen::Vector1d::Constant( 10.0 ), Eigen::Vector1d::Constant( 20.0 ), Eigen::Vector1d::Constant( 30.0 ) },
            { 1.0, 2.0, 3.0 },
            receiver );
    const std::vector< unsigned int > originalIds = dataset.getObservationIdsForSet( setId );
    dataset.removeObservations( ObservationSelectionCondition<>::timeBounds( 1.0, 1.0 ) );

    checkIds( dataset.getObservationIdsForSet( setId ), { originalIds.at( 1 ), originalIds.at( 2 ) } );
    BOOST_CHECK_THROW( dataset.getObservationRow( originalIds.front( ) ), std::out_of_range );
    BOOST_CHECK_NO_THROW( dataset.getObservationValue( originalIds.back( ) ) );

    dataset.addObservationsToSet( setId, { Eigen::Vector1d::Constant( 40.0 ) }, { 0.0 } );
    const std::vector< unsigned int > updatedIds = dataset.getObservationIdsForSet( setId );
    BOOST_CHECK_GT( updatedIds.front( ), originalIds.back( ) );
    checkIds( std::vector< unsigned int >( updatedIds.begin( ) + 1, updatedIds.end( ) ), { originalIds.at( 1 ), originalIds.at( 2 ) } );
}












BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
