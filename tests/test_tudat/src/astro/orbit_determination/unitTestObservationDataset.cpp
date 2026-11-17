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

//! Verify weighted design-matrix output applies sparse observation weights.
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
                                                                                   0.0,
                                                                                   std::vector< Eigen::VectorXd >( ),
                                                                                   sparseWeights );

    const Eigen::MatrixXd lowerWeightSquareRoot = Eigen::LLT< Eigen::MatrixXd >( denseWeights ).matrixL( );
    const Eigen::MatrixXd expectedWeightedDesignMatrix = lowerWeightSquareRoot.transpose( ) * normalizedDesignMatrix;

    BOOST_CHECK( covarianceOutput.hasFullWeightMatrix( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( covarianceOutput.getWeightsMatrix( ).toDense( ), denseWeights, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( covarianceOutput.getNormalizedWeightedDesignMatrix( ), expectedWeightedDesignMatrix, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( covarianceOutput.getUnnormalizedWeightedDesignMatrix( ), expectedWeightedDesignMatrix, 1.0E-15 );
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
//! Verify row conditions select links, values, status, and dependent variables correctly.
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
    checkIds( dataset.getObservationIdsMatchingCondition( ObservationSelectionCondition< double, double >::setId( 1 ) ), { 3, 4 } );

    // Link-based conditions must match rows by complete link definition, link-end type and individual link end.
    checkIds( dataset.getObservationIdsMatchingCondition(
                      ObservationSelectionCondition< double, double >::linkDefinition( station2LinkDefinition ) ),
              { 3, 4 } );
    checkIds( dataset.getObservationIdsMatchingCondition( ObservationSelectionCondition< double, double >::linkEndType( transmitter ) ),
              { 0, 1, 2, 3, 4 } );
    checkIds( dataset.getObservationIdsMatchingCondition( ObservationSelectionCondition< double, double >::linkEndType( reflector1 ) ),
              {} );
    checkIds( dataset.getObservationIdsMatchingCondition(
                      ObservationSelectionCondition< double, double >::linkEnd( transmitter, LinkEndId( "Earth", "Station1" ) ) ),
              { 0, 1, 2 } );

    // One-sided time comparison conditions must select expected rows.
    checkIds( dataset.getObservationIdsMatchingCondition( ObservationSelectionCondition< double, double >::timeGreaterEqual( 4.0 ) ),
              { 3, 4 } );
    checkIds( dataset.getObservationIdsMatchingCondition( ObservationSelectionCondition< double, double >::timeGreaterThan( 4.0 ) ),
              { 4 } );
    checkIds( dataset.getObservationIdsMatchingCondition( ObservationSelectionCondition< double, double >::timeLessEqual( 2.0 ) ),
              { 0, 1 } );
    checkIds( dataset.getObservationIdsMatchingCondition( ObservationSelectionCondition< double, double >::timeLessThan( 2.0 ) ), { 0 } );

    // Boolean condition composition must support unions and complements of row selections.
    checkIds(
            dataset.getObservationIdsMatchingCondition( ObservationSelectionCondition< double, double >::observableType( one_way_range ) ||
                                                        ObservationSelectionCondition< double, double >::timeBounds( 4.5, 5.5 ) ),
            { 0, 1, 2, 4 } );
    checkIds( dataset.getObservationIdsMatchingCondition( !ObservationSelectionCondition< double, double >::timeBounds( 2.0, 4.0 ) ),
              { 0, 4 } );

    // The kind, child conditions and values of a condition can be retrieved.
    const ObservationSelectionCondition< double, double > rangeOrLateCondition =
            ObservationSelectionCondition< double, double >::observableType( one_way_range ) ||
            ObservationSelectionCondition< double, double >::timeGreaterEqual( 4.0 );
    BOOST_CHECK( ( rangeOrLateCondition.getType( ) == ObservationSelectionCondition< double, double >::Type::or_condition ) );
    BOOST_CHECK_EQUAL( rangeOrLateCondition.getChildConditions( ).size( ), 2 );
    BOOST_CHECK( rangeOrLateCondition.getChildConditions( ).at( 0 ).getObservableType( ) == one_way_range );
    BOOST_CHECK( ( rangeOrLateCondition.getChildConditions( ).at( 1 ).getType( ) ==
                   ObservationSelectionCondition< double, double >::Type::time_greater_equal ) );
    BOOST_CHECK_EQUAL( rangeOrLateCondition.getChildConditions( ).at( 1 ).getTimeValue( ), 4.0 );

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

    const Eigen::VectorXd elevationLimit = Eigen::VectorXd::Constant( 1, 0.3 );
    checkIds( dataset.getObservationIdsMatchingCondition(
                      ObservationSelectionCondition< double, double >::dependentVariableGreaterThan( elevationSettings, elevationLimit ) ),
              { 1 } );

    // A less-than condition selects only rows that store the dependent variable, unlike the negated greater-than condition.
    checkIds( dataset.getObservationIdsMatchingCondition(
                      ObservationSelectionCondition< double, double >::dependentVariableLessThan( elevationSettings, elevationLimit ) ),
              { 0, 2 } );
    checkIds( dataset.getObservationIdsMatchingCondition(
                      !ObservationSelectionCondition< double, double >::dependentVariableGreaterThan( elevationSettings, elevationLimit ) ),
              { 0, 2, 3, 4 } );

    // Length-one limits broadcast across vector-valued observations and select a row when any component exceeds the limit.
    checkIds( dataset.getObservationIdsMatchingCondition(
                      ObservationSelectionCondition< double, double >::observableType( angular_position ) &&
                      ObservationSelectionCondition< double, double >::residualAbsoluteValueGreaterThan(
                              Eigen::VectorXd::Constant( 1, 0.6 ) ) ),
              { 4 } );
    checkIds( dataset.getObservationIdsMatchingCondition(
                      ObservationSelectionCondition< double, double >::observableType( angular_position ) &&
                      ObservationSelectionCondition< double, double >::observationAbsoluteValueGreaterThan(
                              Eigen::VectorXd::Constant( 1, 3.5 ) ) ),
              { 4 } );

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

    // With several compatible dependent variables, the first one can be selected explicitly.
    checkIds( ambiguousDependentVariableDataset.getObservationIdsMatchingCondition(
                      ObservationSelectionCondition< double, double >::dependentVariableGreaterThan(
                              elevationSettings, elevationLimit, true ) ),
              { 1 } );

    // A single dependent variable is retrieved by its settings, and is empty for sets that do not store it.
    const std::vector< Eigen::VectorXd > elevationValues = dataset.getDependentVariableValues(
            ObservationSelectionCondition< double, double >::all( ), ObservationOrdering::internal, elevationSettings );
    BOOST_CHECK_EQUAL( elevationValues.size( ), 5 );
    BOOST_CHECK_EQUAL( elevationValues.at( 1 )( 0 ), 0.4 );
    BOOST_CHECK_EQUAL( elevationValues.at( 4 ).size( ), 0 );
    BOOST_CHECK_EQUAL( ambiguousDependentVariableDataset.getSingleDependentVariable( 1, elevationSettings, true )( 0 ), 0.4 );

    ObservationDataset< double, double > noDependentVariableBookkeepingDataset;
    noDependentVariableBookkeepingDataset.addObservationSet(
            one_way_range,
            station1LinkDefinition,
            { ( Eigen::VectorXd( 1 ) << 10.0 ).finished( ), ( Eigen::VectorXd( 1 ) << 20.0 ).finished( ) },
            { 1.0, 2.0 },
            receiver );
    checkIds( noDependentVariableBookkeepingDataset.getObservationIdsMatchingCondition(
                      ObservationSelectionCondition< double, double >::dependentVariableGreaterThan( elevationSettings, elevationLimit ) ),
              {} );

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
 * Verifies that snapshot ordered vector data follows dataset ordered-output order.
 *
 * Test outline: inserts angular-position rows before one-way range rows, so
 * dataset row order differs from ordered observation-vector-data order. Snapshot
 * selections retain dataset order while numerical vector data use the established
 * observable/link/set ordering, including after interleaved appends.
 */
//! Verify ordered vector snapshots reorder selected rows without changing membership.
BOOST_AUTO_TEST_CASE( test_dataset_snapshot_ordered_vector_data_reorders_selected_rows )
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

    const auto snapshot = dataset.getRows( );

    // The snapshot selection itself is in dataset row order: angular rows first, then range rows.
    checkIds( dataset.getObservationIds( ), { 0, 1, 2, 3 } );

    // Estimation projects selected identities into the established numerical order.
    checkIds( dataset.getObservationIds( ObservationSelectionCondition<>::all( ), ObservationOrdering::estimation ), { 2, 3, 0, 1 } );

    // The numerical vector data use the same ordering: range before angular position.
    checkIds( dataset.createObservationVectorData( true ).getObservationIds( ), { 2, 3, 0, 0, 1, 1 } );

    // Appending to the first set after another set exists makes storage noncontiguous.
    dataset.addObservationsToSet( 0, { Eigen::Vector2d( 24.0, 25.0 ) }, { 2.5 }, {}, {}, {}, true );
    BOOST_CHECK_EQUAL( snapshot.size( ), 4 );
    BOOST_CHECK_EQUAL( snapshot.at( 0 ).observationId_, 0 );
    const auto appendedVectorData = dataset.createObservationVectorData( );
    checkIds( appendedVectorData.getObservationIds( ), { 2, 3, 4, 4, 0, 0, 1, 1 } );
    checkIds( appendedVectorData.getScalarComponentIds( ), { 4, 5, 6, 7, 0, 1, 2, 3 } );
    Eigen::VectorXd expectedAppended( 8 );
    expectedAppended << 10.0, 11.0, 24.0, 25.0, 20.0, 21.0, 22.0, 23.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( appendedVectorData.getObservationVector( ), expectedAppended, 1.0E-15 );
    BOOST_CHECK_EQUAL( appendedVectorData.getVectorRow( 4, 1 ), 3 );
    BOOST_CHECK_EQUAL( appendedVectorData.getVectorRow( 0, 1 ), 5 );
}

/*!
 * Verifies snapshots, rejection/restoration and reduced dataset creation.
 *
 * Test outline: builds a dataset, selects rows with conditions, creates snapshots and
 * reduced datasets, rejects/restores observations and checks observation vector data sizes.
 * It also confirms that snapshots survive structural mutations.
 */
//! Verify condition snapshots drive rejection and reduced-dataset creation.
BOOST_AUTO_TEST_CASE( test_dataset_condition_snapshot_rejection_and_reduced_dataset )
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

    const auto snapshot = dataset.getObservations( middleTimes );

    // A snapshot must project only rows that satisfy its condition.
    BOOST_CHECK_EQUAL( snapshot.size( ), 2 );
    Eigen::VectorXd expectedSnapshotObservations( 2 );
    expectedSnapshotObservations << 20.0, 30.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( concatenateInspectionValues( snapshot ), expectedSnapshotObservations, 1.0E-15 );

    const std::shared_ptr< ObservationDataset< double, double > > keptDataset = dataset.createNewAndKeep( middleTimes );

    // A reduced keep-dataset must contain only selected rows and preserve their scalar observation vector data.
    BOOST_CHECK_EQUAL( keptDataset->getNumberOfObservationSets( ), 1 );
    BOOST_CHECK_EQUAL( keptDataset->getNumberOfObservations( ), 2 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            keptDataset->createObservationVectorData( ).getObservationVector( ), expectedSnapshotObservations, 1.0E-15 );

    const std::shared_ptr< ObservationDataset< double, double > > droppedDataset = dataset.createNewAndDrop( middleTimes );
    Eigen::VectorXd expectedDroppedObservations( 2 );
    expectedDroppedObservations << 10.0, 40.0;

    // A reduced drop-dataset must contain the complement of the selected rows.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            droppedDataset->createObservationVectorData( ).getObservationVector( ), expectedDroppedObservations, 1.0E-15 );

    dataset.rejectObservations( ObservationSelectionCondition< double, double >::timeBounds( 2.5, 3.5 ), "test rejection" );

    // Rejection must store row status and reason while removing the row from active observation vector data.
    BOOST_CHECK( !dataset.getObservationRow( 2 ).isActive_ );
    BOOST_CHECK_EQUAL( dataset.getObservationRow( 2 ).rejectionReason_, "test rejection" );

    Eigen::VectorXd expectedActiveObservations( 3 );
    expectedActiveObservations << 10.0, 20.0, 40.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            dataset.createObservationVectorData( ).getObservationVector( ), expectedActiveObservations, 1.0E-15 );
    BOOST_CHECK_EQUAL( dataset.createObservationVectorData( true ).getObservationVector( ).size( ), 4 );
    BOOST_CHECK_EQUAL( dataset.getObservationIds( ObservationSelectionCondition< double, double >::rejected( ) ).size( ), 1 );

    dataset.restoreObservations( ObservationSelectionCondition< double, double >::rejected( ) );

    // Restoration must reactivate rejected rows and clear their rejection reason.
    BOOST_CHECK( dataset.getObservationRow( 2 ).isActive_ );
    BOOST_CHECK( dataset.getObservationRow( 2 ).rejectionReason_.empty( ) );
    BOOST_CHECK_EQUAL( dataset.createObservationVectorData( ).getObservationVector( ).size( ), 4 );

    dataset.rejectObservations( ObservationSelectionCondition< double, double >::timeBounds( 2.5, 3.5 ), "delete rejection" );
    dataset.deleteRejectedObservations( );

    // Deleting rejected rows is explicit and physically removes rows from the dataset.
    Eigen::VectorXd expectedAfterRejectedDeletion( 3 );
    expectedAfterRejectedDeletion << 10.0, 20.0, 40.0;
    BOOST_CHECK_EQUAL( dataset.getNumberOfObservations( ), 3 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            dataset.createObservationVectorData( true ).getObservationVector( ), expectedAfterRejectedDeletion, 1.0E-15 );
    BOOST_CHECK_EQUAL( dataset.getObservationIds( ObservationSelectionCondition< double, double >::rejected( ) ).size( ), 0 );

    dataset.removeObservations( ObservationSelectionCondition< double, double >::timeBounds( 3.5, 4.5 ) );

    // Direct condition-based removal must provide the destructive counterpart to createNewAndDrop.
    Eigen::VectorXd expectedAfterConditionRemoval( 2 );
    expectedAfterConditionRemoval << 10.0, 20.0;
    BOOST_CHECK_EQUAL( dataset.getNumberOfObservations( ), 2 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            dataset.createObservationVectorData( true ).getObservationVector( ), expectedAfterConditionRemoval, 1.0E-15 );

    const auto invalidatedSnapshot = dataset.getRows( );
    dataset.addObservationSet( one_way_range, linkDefinition, std::vector< Eigen::VectorXd >( ), std::vector< double >( ), receiver );

    // Structural dataset mutations leave previously extracted snapshots intact.
    BOOST_CHECK_EQUAL( invalidatedSnapshot.size( ), 2 );

    ObservationDataset< double, double > linkMutationDataset;
    linkMutationDataset.addObservationSet( one_way_range, linkDefinition, { Eigen::Vector1d::Constant( 1.0 ) }, { 1.0 }, receiver );
    const auto linkMutationSnapshot = linkMutationDataset.getRows( );
    linkMutationDataset.setLinkEndReferencePoint( "Earth", "StationX", transmitter );

    // Updating link-end reference points leaves snapshot membership intact.
    BOOST_CHECK_EQUAL( linkMutationSnapshot.size( ), 1 );
    BOOST_CHECK_EQUAL( linkMutationDataset.getLinkDefinition( linkMutationDataset.getObservationSetMetadata( 0 ).linkDefinitionId_ )
                               .linkEnds_.at( transmitter )
                               .getReferencePointName( ),
                       "StationX" );
}

//! Verify snapshots outlive datasets and dataset assignment.
BOOST_AUTO_TEST_CASE( test_dataset_snapshot_lifetime )
{
    const LinkDefinition station1LinkDefinition = createOneWayLinkDefinition( "Station1" );

    const auto expiredSnapshot = [ station1LinkDefinition ]( ) {
        ObservationDataset< double, double > localDataset;
        localDataset.addObservationSet( one_way_range, station1LinkDefinition, { Eigen::Vector1d::Constant( 1.0 ) }, { 1.0 }, receiver );
        return localDataset.getObservations( );
    }( );
    BOOST_CHECK_EQUAL( expiredSnapshot.size( ), 1 );

    ObservationDataset< double, double > assignedDataset;
    assignedDataset.addObservationSet( one_way_range, station1LinkDefinition, { Eigen::Vector1d::Constant( 1.0 ) }, { 1.0 }, receiver );
    const auto preAssignmentSnapshot = assignedDataset.getRows( );
    ObservationDataset< double, double > replacementDataset;
    replacementDataset.addObservationSet( one_way_range, station1LinkDefinition, { Eigen::Vector1d::Constant( 2.0 ) }, { 2.0 }, receiver );
    assignedDataset = replacementDataset;
    BOOST_CHECK_EQUAL( preAssignmentSnapshot.size( ), 1 );
}

/*!
 * Verifies compact, block and sparse off-diagonal weight handling.
 *
 * Test outline: covers scalar weights, per-observation blocks, set-level blocks,
 * observation-id selected cross blocks, symmetric transpose insertion,
 * component selection and rejection/restoration of weighted observations. It
 * checks both compact weight vectors and materialized sparse observation vector data
 * matrices.
 */
//! Verify compact weight policies and full matrix weights produce the expected coefficients.
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
            constantScalarWeights( 5.0 ) );

    // Scalar-weight construction must remain diagonal-only and expose no per-observation matrix blocks.
    BOOST_CHECK_EQUAL( scalarSetId, 0 );
    Eigen::VectorXd expectedScalarWeights( 4 );
    expectedScalarWeights << 5.0, 5.0, 5.0, 5.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            scalarWeightDataset.getWeightDiagonal( ObservationSelectionCondition< double, double >::setId( scalarSetId ) ),
            expectedScalarWeights,
            1.0E-15 );
    const Eigen::MatrixXd expectedScalarWeightMatrix = 5.0 * Eigen::MatrixXd::Identity( 4, 4 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            scalarWeightDataset.createObservationVectorData( ).getSparseWeightMatrix( ).toDense( ), expectedScalarWeightMatrix, 1.0E-15 );
    BOOST_CHECK_EQUAL( scalarWeightDataset.createObservationVectorData( ).getSparseWeightMatrix( ).nonZeros( ), 4 );
    BOOST_CHECK( !scalarWeightDataset.createObservationVectorData( ).hasOffDiagonalWeights( ) );
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
            constantBlockWeights( observationWeightBlock ) );
    Eigen::MatrixXd expectedObservationBlockMatrix = Eigen::MatrixXd::Zero( 4, 4 );
    expectedObservationBlockMatrix.block( 0, 0, 2, 2 ) = observationWeightBlock;
    expectedObservationBlockMatrix.block( 2, 2, 2, 2 ) = observationWeightBlock;

    // Per-observation weight blocks must be repeated for each vector-valued observation and force sparse weighting.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( blockWeightDataset.createObservationVectorData( ).getSparseWeightMatrix( ).toDense( ),
                                       expectedObservationBlockMatrix,
                                       1.0E-15 );
    BOOST_CHECK_EQUAL( blockWeightDataset.createObservationVectorData( ).getSparseWeightMatrix( ).nonZeros( ), 8 );
    BOOST_CHECK( blockWeightDataset.createObservationVectorData( ).hasOffDiagonalWeights( ) );
    BOOST_CHECK( blockWeightDataset.hasWeightMatrixForObservation( blockWeightDataset.getObservationIdsForSet( blockSetId ).at( 0 ) ) );
    Eigen::VectorXd expectedBlockWeightDiagonal( 4 );
    expectedBlockWeightDiagonal << 2.0, 3.0, 2.0, 3.0;

    // Compact weight-vector access must return the diagonal of repeated per-observation blocks.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            blockWeightDataset.getWeightDiagonal( ObservationSelectionCondition< double, double >::setId( blockSetId ) ),
            expectedBlockWeightDiagonal,
            1.0E-15 );

    ObservationDataset< double, double > setConstantMatrixDataset;
    const int constantMatrixSetId = setConstantMatrixDataset.addObservationSet(
            angular_position,
            linkDefinition,
            { ( Eigen::Vector2d( ) << 10.0, 11.0 ).finished( ), ( Eigen::Vector2d( ) << 20.0, 21.0 ).finished( ) },
            { 1.0, 2.0 },
            receiver );
    setConstantMatrixDataset.setConstantSingleObservationMatrixWeight(
            ObservationSelectionCondition< double, double >::setId( constantMatrixSetId ), observationWeightBlock );
    BOOST_CHECK( setConstantMatrixDataset.getObservationSetMetadata( constantMatrixSetId ).weightStructure_ ==
                 ObservationWeightStructure::per_observation );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( setConstantMatrixDataset.createObservationVectorData( ).getSparseWeightMatrix( ).toDense( ),
                                       expectedObservationBlockMatrix,
                                       1.0E-15 );

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
            conditionWeightDataset.getWeightDiagonal( ObservationSelectionCondition< double, double >::setId( conditionRangeSetId ) ),
            ( 7.0 * Eigen::Vector2d::Ones( ) ),
            1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( conditionWeightDataset.getWeightMatrixForObservation(
                                               conditionWeightDataset.getObservationIdsForSet( conditionAngularSetId ).at( 0 ) ),
                                       observationWeightBlock,
                                       1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            conditionWeightDataset.getWeightValue( conditionWeightDataset.getObservationIdsForSet( conditionAngularSetId ).at( 1 ) ),
            ( Eigen::Vector2d( ) << 8.0, 9.0 ).finished( ),
            1.0E-15 );

    Eigen::MatrixXd setWeightBlock = Eigen::MatrixXd::Zero( 4, 4 );
    setWeightBlock << 1.0, 0.1, 0.2, 0.3, 0.1, 2.0, 0.4, 0.5, 0.2, 0.4, 3.0, 0.6, 0.3, 0.5, 0.6, 4.0;
    ObservationDataset< double, double > setBlockWeightDataset;
    const int setBlockSetId = setBlockWeightDataset.addObservationSetWithWeights(
            angular_position,
            linkDefinition,
            { ( Eigen::Vector2d( ) << 10.0, 11.0 ).finished( ), ( Eigen::Vector2d( ) << 20.0, 21.0 ).finished( ) },
            { 1.0, 2.0 },
            receiver,
            setBlockWeights( setWeightBlock ) );

    // A set-level weight block must be stored as one matrix spanning the full scalar range of the set.
    BOOST_CHECK( setBlockWeightDataset.hasWeightMatrixForSet( setBlockSetId ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( setBlockWeightDataset.getWeightMatrixForSet( setBlockSetId ), setWeightBlock, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            setBlockWeightDataset.createObservationVectorData( ).getSparseWeightMatrix( ).toDense( ), setWeightBlock, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            setBlockWeightDataset.createObservationVectorData( ).getWeightVector( ), setWeightBlock.diagonal( ), 1.0E-15 );

    ObservationDataset< double, double > extraBlockDataset;
    extraBlockDataset.addObservationSetWithWeights(
            one_way_range,
            linkDefinition,
            { Eigen::Vector1d::Constant( 10.0 ), Eigen::Vector1d::Constant( 20.0 ), Eigen::Vector1d::Constant( 30.0 ) },
            { 1.0, 2.0, 3.0 },
            receiver,
            constantScalarWeights( 1.0 ) );
    const std::vector< unsigned int > extraBlockObservationIds = extraBlockDataset.getObservationIdsForSet( 0 );
    const Eigen::Matrix2d extraWeightBlock = ( Eigen::Matrix2d( ) << 4.0, 0.25, 0.25, 5.0 ).finished( );
    extraBlockDataset.setWeightBlock( { extraBlockObservationIds.at( 0 ), extraBlockObservationIds.at( 2 ) },
                                      { extraBlockObservationIds.at( 0 ), extraBlockObservationIds.at( 2 ) },
                                      extraWeightBlock,
                                      {},
                                      {} );

    // Arbitrary observation-id blocks update the same effective symmetric matrix.
    BOOST_CHECK( extraBlockDataset.hasExtraWeightBlocks( ) );
    Eigen::MatrixXd expectedExtraBlockWeightMatrix = Eigen::MatrixXd::Identity( 3, 3 );
    expectedExtraBlockWeightMatrix( 0, 0 ) = 4.0;
    expectedExtraBlockWeightMatrix( 2, 2 ) = 5.0;
    expectedExtraBlockWeightMatrix( 0, 2 ) = 0.25;
    expectedExtraBlockWeightMatrix( 2, 0 ) = 0.25;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            extraBlockDataset.createObservationVectorData( ).getSparseWeightMatrix( ).toDense( ), expectedExtraBlockWeightMatrix, 1.0E-15 );

    ObservationDataset< double, double > symmetricComponentBlockDataset;
    symmetricComponentBlockDataset.addObservationSetWithWeights( angular_position,
                                                                 linkDefinition,
                                                                 { ( Eigen::Vector2d( ) << 10.0, 11.0 ).finished( ),
                                                                   ( Eigen::Vector2d( ) << 20.0, 21.0 ).finished( ),
                                                                   ( Eigen::Vector2d( ) << 30.0, 31.0 ).finished( ) },
                                                                 { 1.0, 2.0, 3.0 },
                                                                 receiver,
                                                                 constantScalarWeights( 1.0 ) );
    const std::vector< unsigned int > symmetricObservationIds = symmetricComponentBlockDataset.getObservationIdsForSet( 0 );
    const Eigen::MatrixXd crossComponentBlock = ( Eigen::Matrix< double, 2, 1 >( ) << 0.7, 0.8 ).finished( );
    symmetricComponentBlockDataset.setWeightBlock( { symmetricObservationIds.at( 0 ), symmetricObservationIds.at( 1 ) },
                                                   { symmetricObservationIds.at( 2 ) },
                                                   crossComponentBlock,
                                                   { 0 },
                                                   { 1 } );

    // Component-selected blocks must populate both symmetric matrix entries.
    Eigen::MatrixXd expectedSymmetricComponentWeightMatrix = Eigen::MatrixXd::Identity( 6, 6 );
    expectedSymmetricComponentWeightMatrix( 0, 5 ) = 0.7;
    expectedSymmetricComponentWeightMatrix( 2, 5 ) = 0.8;
    expectedSymmetricComponentWeightMatrix( 5, 0 ) = 0.7;
    expectedSymmetricComponentWeightMatrix( 5, 2 ) = 0.8;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( symmetricComponentBlockDataset.createObservationVectorData( ).getSparseWeightMatrix( ).toDense( ),
                                       expectedSymmetricComponentWeightMatrix,
                                       1.0E-15 );

    ObservationDataset< double, double > rejectedSetBlockDataset;
    Eigen::MatrixXd fullSetWeightBlock = Eigen::MatrixXd::Zero( 6, 6 );
    for( int row = 0; row < fullSetWeightBlock.rows( ); ++row )
    {
        for( int column = 0; column < fullSetWeightBlock.cols( ); ++column )
        {
            fullSetWeightBlock( row, column ) = row == column
                    ? static_cast< double >( row + 1 )
                    : 0.01 * static_cast< double >( 10 * ( std::min( row, column ) + 1 ) + std::max( row, column ) + 1 );
        }
    }
    rejectedSetBlockDataset.addObservationSetWithWeights( angular_position,
                                                          linkDefinition,
                                                          { ( Eigen::Vector2d( ) << 10.0, 11.0 ).finished( ),
                                                            ( Eigen::Vector2d( ) << 20.0, 21.0 ).finished( ),
                                                            ( Eigen::Vector2d( ) << 30.0, 31.0 ).finished( ) },
                                                          { 1.0, 2.0, 3.0 },
                                                          receiver,
                                                          setBlockWeights( fullSetWeightBlock ) );
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
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( rejectedSetBlockDataset.createObservationVectorData( ).getSparseWeightMatrix( ).toDense( ),
                                       expectedRejectedWeightBlock,
                                       1.0E-15 );
    rejectedSetBlockDataset.restoreObservations( ObservationSelectionCondition< double, double >::rejected( ) );

    // Restoring rejected rows must restore the full original weight block in the observation vector data.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            rejectedSetBlockDataset.createObservationVectorData( ).getSparseWeightMatrix( ).toDense( ), fullSetWeightBlock, 1.0E-15 );
}

/*!
 * Verifies the precedence between set-level, per-observation and extra block weights.
 *
 * Test outline: creates one vector-valued set with a full set-level block, then
 * overwrites part of it with an explicit per-observation block and finally
 * overwrites cross-observation entries with an extra scalar-component block. It
 * checks the resulting dense sparse matrix and confirms that conflicts report
 * the affected observation-vector and scalar-component indices.
 */
//! Verify weight assignments replace only the addressed matrix entries.
BOOST_AUTO_TEST_CASE( test_dataset_weight_assignments_replace_addressed_entries )
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
    perObservationWeightBlock << 10.0, 0.7, 0.7, 11.0;
    dataset.setWeightMatrixForObservation( observationIds.at( 0 ), perObservationWeightBlock );
    const Eigen::Vector2d perObservationDiagonalWeights = ( Eigen::Vector2d( ) << 12.0, 13.0 ).finished( );
    dataset.setConstantSingleObservationDiagonalWeight( ObservationSelectionCondition< double, double >::timeBounds( 1.5, 2.5 ),
                                                        perObservationDiagonalWeights );
    BOOST_CHECK( dataset.getObservationSetMetadata( setId ).weightStructure_ == ObservationWeightStructure::per_set );

    Eigen::Matrix2d extraWeightBlock;
    extraWeightBlock << 20.0, 21.0, 22.0, 23.0;
    dataset.setWeightBlock( { observationIds.at( 0 ) }, { observationIds.at( 1 ) }, extraWeightBlock );

    const auto observationVectorData = dataset.createObservationVectorData( );

    Eigen::MatrixXd expectedWeightMatrix = setWeightBlock;
    expectedWeightMatrix.block( 0, 0, 2, 2 ) = perObservationWeightBlock;
    expectedWeightMatrix.block( 2, 2, 2, 2 ) = perObservationDiagonalWeights.asDiagonal( );
    expectedWeightMatrix.block( 0, 2, 2, 2 ) = extraWeightBlock;
    expectedWeightMatrix.block( 2, 0, 2, 2 ) = extraWeightBlock.transpose( );

    // Set-level blocks provide the baseline, explicit matrix/vector observation weights overwrite that baseline,
    // and extra scalar-component blocks overwrite both lower-priority layers while preserving matrix symmetry.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( observationVectorData.getSparseWeightMatrix( ).toDense( ), expectedWeightMatrix, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( observationVectorData.getWeightVector( ), expectedWeightMatrix.diagonal( ), 1.0E-15 );

    // A later set-wide assignment replaces previous row and cross-row coefficients as well.
    dataset.setWeightMatrixForSet( setId, setWeightBlock );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            dataset.createObservationVectorData( ).getSparseWeightMatrix( ).toDense( ), setWeightBlock, 1.0E-15 );
}

/*!
 * Verifies that later assignments consistently replace appended weights.
 *
 * Test outline: appends rows before setting a full set-level weight block. The
 * first dataset omits append weights, and the second supplies them explicitly.
 * A later set-level assignment must replace the addressed matrix in both cases.
 */
//! Verify appended default weights do not overwrite a subsequently assigned set block.
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

    // The later set assignment replaces appended defaults.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            implicitAppendDataset.createObservationVectorData( ).getSparseWeightMatrix( ).toDense( ), setWeightBlock, 1.0E-15 );

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

    // The later full-set assignment also replaces explicitly supplied append weights.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            explicitAppendDataset.createObservationVectorData( ).getSparseWeightMatrix( ).toDense( ), setWeightBlock, 1.0E-15 );
}

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
 * sparse weights and compares the returned parameter update and normal matrix
 * against an independently assembled dense reference calculation.
 */
//! Verify least-squares adjustment uses a sparse full observation weight matrix.
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
            linear_algebra::performLeastSquaresAdjustmentFromDesignMatrix(
                    designMatrix, residuals, sparseWeights, Eigen::MatrixXd::Zero( designMatrix.cols( ), designMatrix.cols( ) ) );

    // Sparse-weight least squares must produce the same normal matrix and parameter update as the dense reference.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( sparseLeastSquaresOutput.second, expectedNormalMatrix, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( sparseLeastSquaresOutput.first, expectedParameterUpdate, 1.0E-15 );
}

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

//! Verify residual writeback rejects vector data created before the dataset was changed.
BOOST_AUTO_TEST_CASE( test_observation_vector_data_writeback_requires_unchanged_dataset )
{
    ObservationDataset<> dataset;
    const unsigned int setId = dataset.addObservationSet( one_way_range,
                                                          createOneWayLinkDefinition( "Station1" ),
                                                          { Eigen::Vector1d::Constant( 10.0 ), Eigen::Vector1d::Constant( 20.0 ) },
                                                          { 1.0, 2.0 },
                                                          receiver );
    const auto vectorData = dataset.createObservationVectorData( );
    const Eigen::Vector2d residuals( 3.0, 4.0 );

    // A metadata change can alter design-matrix ordering without changing vector dimensions.
    dataset.setLinkEndReferencePoint( "Earth", "Station2", transmitter );
    BOOST_CHECK_THROW( dataset.setResidualVector( vectorData, residuals ), std::runtime_error );

    const auto currentVectorData = dataset.createObservationVectorData( );
    dataset.setResidualVector( currentVectorData, residuals );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( dataset.getResidualVectorForSet( setId ), residuals, 1.0E-15 );
}

//! Preserve the supplied block structure even when the numerical matrices happen to be diagonal.
BOOST_AUTO_TEST_CASE( test_weight_structure_from_settings )
{
    using Structure = ObservationWeightStructure;
    const std::vector< std::pair< std::shared_ptr< ObservationWeightSettings >, Structure > > settings = {
        { constantScalarWeights( 2.0 ), Structure::diagonal },
        { scalarPerObservationWeights( { 2.0, 3.0 } ), Structure::diagonal },
        { diagonalPerObservationWeights( { Eigen::Vector2d::Ones( ), Eigen::Vector2d::Ones( ) } ), Structure::diagonal },
        { constantBlockWeights( Eigen::Matrix2d::Identity( ) ), Structure::per_observation },
        { blockPerObservationWeights( { Eigen::Matrix2d::Identity( ), Eigen::Matrix2d::Identity( ) } ), Structure::per_observation },
        { setBlockWeights( Eigen::Matrix4d::Identity( ) ), Structure::per_set }
    };
    for( const auto& entry : settings )
    {
        ObservationDataset<> dataset;
        const unsigned int setId = dataset.addObservationSetWithWeights( angular_position,
                                                                         createOneWayLinkDefinition( "Station1" ),
                                                                         { Eigen::Vector2d::Ones( ), Eigen::Vector2d::Ones( ) },
                                                                         { 1.0, 2.0 },
                                                                         receiver,
                                                                         entry.first );
        BOOST_CHECK( dataset.getObservationSetMetadata( setId ).weightStructure_ == entry.second );
        BOOST_CHECK( dataset.createObservationVectorData( ).getSetMetadata( setId ).weightStructure_ == entry.second );
        const ObservationDataset<> copied( dataset );
        BOOST_CHECK( copied.getObservationSetMetadata( setId ).weightStructure_ == entry.second );
        ObservationDataset<> imported;
        const unsigned int importedSet = imported.addObservationSetFromDataset( dataset, setId );
        BOOST_CHECK( imported.getObservationSetMetadata( importedSet ).weightStructure_ == entry.second );
    }
}

//! Update both sets for cross-set entries without confusing zero entries or diagonal replacements with correlations.
BOOST_AUTO_TEST_CASE( test_weight_structure_after_weight_changes )
{
    using Structure = ObservationWeightStructure;
    ObservationDataset<> dataset;
    const auto link = createOneWayLinkDefinition( "Station1" );
    const unsigned int firstSet = dataset.addObservationSet(
            angular_position, link, { Eigen::Vector2d::Ones( ), Eigen::Vector2d::Ones( ) }, { 1.0, 2.0 }, receiver );
    const unsigned int secondSet = dataset.addObservationSet( angular_position, link, { Eigen::Vector2d::Ones( ) }, { 3.0 }, receiver );
    const auto firstIds = dataset.getObservationIdsForSet( firstSet );
    const auto secondIds = dataset.getObservationIdsForSet( secondSet );
    const Eigen::MatrixXd zero = Eigen::MatrixXd::Zero( 1, 1 );
    const Eigen::MatrixXd correlation = Eigen::MatrixXd::Constant( 1, 1, 0.1 );
    dataset.setWeightBlock( { firstIds.front( ) }, secondIds, zero, { 1 }, { 0 } );
    BOOST_CHECK( dataset.getObservationSetMetadata( firstSet ).weightStructure_ == Structure::diagonal );
    BOOST_CHECK( dataset.getObservationSetMetadata( secondSet ).weightStructure_ == Structure::diagonal );
    dataset.setWeightBlock( { firstIds.front( ) }, { firstIds.front( ) }, correlation, { 1 }, { 0 } );
    BOOST_CHECK( dataset.getObservationSetMetadata( firstSet ).weightStructure_ == Structure::per_observation );
    dataset.setWeightBlock( { firstIds.front( ) }, { firstIds.back( ) }, correlation, { 0 }, { 1 } );
    BOOST_CHECK( dataset.getObservationSetMetadata( firstSet ).weightStructure_ == Structure::per_set );

    dataset.setWeightBlock( { firstIds.front( ) }, secondIds, correlation, { 1 }, { 0 } );
    const auto vectorData = dataset.createObservationVectorData( );
    BOOST_CHECK( vectorData.getSetMetadata( firstSet ).weightStructure_ == Structure::inter_set_weights );
    BOOST_CHECK( vectorData.getSetMetadata( secondSet ).weightStructure_ == Structure::inter_set_weights );
    // These assignments do not touch the connection between the two sets.
    dataset.setWeightVectorForSet( firstSet, Eigen::Vector4d::Ones( ) );
    dataset.setWeightMatrixForObservation( secondIds.front( ), Eigen::Matrix2d::Identity( ) );
    BOOST_CHECK( dataset.getObservationSetMetadata( firstSet ).weightStructure_ == Structure::inter_set_weights );
    BOOST_CHECK( dataset.getObservationSetMetadata( secondSet ).weightStructure_ == Structure::inter_set_weights );
    const ObservationDataset<> copied( dataset );
    BOOST_CHECK( copied.getObservationSetMetadata( firstSet ).weightStructure_ == Structure::inter_set_weights );
    ObservationDataset<> imported;
    const unsigned int importedSet = imported.addObservationSetFromDataset( dataset, firstSet );
    BOOST_CHECK( imported.getObservationSetMetadata( importedSet ).weightStructure_ == Structure::per_set );

    dataset.setWeightBlock( { firstIds.front( ) }, secondIds, zero, { 1 }, { 0 } );
    BOOST_CHECK( dataset.getObservationSetMetadata( firstSet ).weightStructure_ == Structure::per_set );
    BOOST_CHECK( dataset.getObservationSetMetadata( secondSet ).weightStructure_ == Structure::per_set );
    // Existing snapshots retain their own metadata.
    BOOST_CHECK( vectorData.getSetMetadata( firstSet ).weightStructure_ == Structure::inter_set_weights );
}

//! Measurement covariance uses complete weights even when the vectors omit rejected observations.
BOOST_AUTO_TEST_CASE( test_complete_observation_covariances )
{
    Eigen::Matrix4d setWeights;
    setWeights << 4.0, 0.0, 1.0, 0.0, 0.0, 9.0, 0.0, 3.0, 1.0, 0.0, 2.0, 0.0, 0.0, 3.0, 0.0, 6.0;
    const Eigen::Matrix2d observationWeights = ( Eigen::Matrix2d( ) << 4.0, 1.0, 1.0, 2.0 ).finished( );
    const Eigen::Matrix2d observationCovariance = ( Eigen::Matrix2d( ) << 2.0, -1.0, -1.0, 4.0 ).finished( ) / 7.0;
    const std::vector< std::shared_ptr< ObservationWeightSettings > > settings = {
        diagonalPerObservationWeights( { Eigen::Vector2d( 4.0, 9.0 ), Eigen::Vector2d( 2.0, 6.0 ) } ),
        constantBlockWeights( observationWeights ),
        setBlockWeights( setWeights )
    };
    const std::vector< Eigen::Matrix2d > firstCovariances = { Eigen::Vector2d( 1.0 / 4.0, 1.0 / 9.0 ).asDiagonal( ),
                                                              observationCovariance,
                                                              Eigen::Vector2d( 2.0 / 7.0, 2.0 / 15.0 ).asDiagonal( ) };
    const std::vector< Eigen::Matrix2d > secondCovariances = { Eigen::Vector2d( 1.0 / 2.0, 1.0 / 6.0 ).asDiagonal( ),
                                                               observationCovariance,
                                                               Eigen::Vector2d( 4.0 / 7.0, 1.0 / 5.0 ).asDiagonal( ) };
    for( std::size_t i = 0; i < settings.size( ); ++i )
    {
        ObservationDataset<> dataset;
        const unsigned int setId = dataset.addObservationSetWithWeights( angular_position,
                                                                         createOneWayLinkDefinition( "Station1" ),
                                                                         { Eigen::Vector2d::Ones( ), Eigen::Vector2d::Ones( ) },
                                                                         { 1.0, 2.0 },
                                                                         receiver,
                                                                         settings.at( i ) );
        const auto ids = dataset.getObservationIdsForSet( setId );
        dataset.rejectObservations( ObservationSelectionCondition<>::timeBounds( 2.0, 2.0 ) );
        const auto selected = dataset.createObservationVectorData( );
        BOOST_CHECK_EQUAL( selected.getObservationVector( ).size( ), 2 );
        dataset.rejectObservations( ObservationSelectionCondition<>::all( ) );
        const auto allRejected = dataset.createObservationVectorData( );
        BOOST_CHECK_EQUAL( allRejected.getObservationVector( ).size( ), 0 );
        // Change the dataset before either snapshot has computed any covariance.
        dataset.setWeightVectorForSet( setId, Eigen::Vector4d::Constant( 8.0 ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                selected.getInverseWeightMatrixForObservation( ids.front( ) ), firstCovariances.at( i ), 1.0E-14 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                selected.getInverseWeightMatrixForObservation( ids.back( ) ), secondCovariances.at( i ), 1.0E-14 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                allRejected.getInverseWeightMatrixForObservation( ids.back( ) ), secondCovariances.at( i ), 1.0E-14 );
        // The ordinary estimation weights still contain only the retained observation's weight block.
        const Eigen::Matrix2d selectedWeights = i == 1 ? observationWeights : Eigen::Vector2d( 4.0, 9.0 ).asDiagonal( ).toDenseMatrix( );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( selected.getSparseWeightMatrix( ).toDense( ), selectedWeights, 1.0E-15 );
    }
}

//! Verify all weight getters expose the same effective matrix coefficients.
BOOST_AUTO_TEST_CASE( test_weight_getters_return_the_effective_matrix )
{
    ObservationDataset<> dataset;
    const unsigned int setId = dataset.addObservationSet( angular_position,
                                                          createOneWayLinkDefinition( "Station1" ),
                                                          { Eigen::Vector2d( 1.0, 2.0 ), Eigen::Vector2d( 3.0, 4.0 ) },
                                                          { 1.0, 2.0 },
                                                          receiver );
    const auto ids = dataset.getObservationIdsForSet( setId );
    const Eigen::Matrix4d setMatrix =
            ( Eigen::Matrix4d( ) << 4.0, 0.1, 0.2, 0.3, 0.1, 5.0, 0.4, 0.5, 0.2, 0.4, 6.0, 0.6, 0.3, 0.5, 0.6, 7.0 ).finished( );
    dataset.setWeightMatrixForSet( setId, setMatrix );
    const Eigen::Matrix2d observationBlock = ( Eigen::Matrix2d( ) << 8.0, 0.7, 0.7, 9.0 ).finished( );
    dataset.setWeightMatrixForObservation( ids.front( ), observationBlock );
    dataset.setWeightBlock( { ids.back( ) }, { ids.back( ) }, Eigen::MatrixXd::Constant( 1, 1, 10.0 ), { 1 }, { 1 } );

    Eigen::Matrix4d expected = setMatrix;
    expected.topLeftCorner< 2, 2 >( ) = observationBlock;
    expected( 3, 3 ) = 10.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            dataset.getWeightDiagonal( ObservationSelectionCondition<>::setId( setId ) ), expected.diagonal( ), 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( dataset.getWeightValue( ids.back( ) ), expected.diagonal( ).tail( 2 ), 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( dataset.getWeightMatrixForSet( setId ), expected, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( dataset.createObservationVectorData( ).getSparseWeightMatrix( ).toDense( ), expected, 1.0E-15 );
}

//! Verify vector data are value snapshots and that residuals can be written back repeatedly within one iteration.
BOOST_AUTO_TEST_CASE( test_observation_vector_data_snapshots_and_repeated_residual_writeback )
{
    ObservationDataset<> dataset;
    auto ancillary = std::make_shared< ObservationAncillarySimulationSettings >( );
    ancillary->setAncillaryDoubleData( doppler_integration_time, 10.0 );
    const auto setId = dataset.addObservationSet( one_way_range,
                                                  createOneWayLinkDefinition( "Station1" ),
                                                  { Eigen::Vector1d::Constant( 10.0 ) },
                                                  { 1.0 },
                                                  receiver,
                                                  {},
                                                  nullptr,
                                                  ancillary );
    const auto snapshot = dataset.getRows( );
    const auto beforeValues = dataset.createObservationVectorData( );
    dataset.setObservationsForSet( setId, { Eigen::Vector1d::Constant( 20.0 ) } );
    BOOST_CHECK_EQUAL( snapshot.size( ), 1 );
    BOOST_CHECK_EQUAL( beforeValues.getObservationVector( )( 0 ), 10.0 );

    const auto beforeAncillary = dataset.createObservationVectorData( );
    dataset.getAncillarySettingsForSet( setId )->setAncillaryDoubleData( doppler_integration_time, 20.0 );
    BOOST_CHECK_EQUAL( beforeAncillary.getAncillarySettingsForSet( setId )->getAncillaryDoubleData( doppler_integration_time ), 10.0 );

    // Residual updates and no-op removals do not invalidate an iteration's mapping.
    const auto current = dataset.createObservationVectorData( );
    dataset.removeObservations( !ObservationSelectionCondition<>::all( ) );
    dataset.setResidualVector( current, Eigen::Vector1d::Constant( 2.0 ) );
    dataset.setResidualVector( current, Eigen::Vector1d::Constant( 3.0 ) );
    BOOST_CHECK_EQUAL( dataset.getResidualValue( 0 )( 0 ), 3.0 );
}

//! Verify design matrices use active mixed-dimension observations in estimation order.
BOOST_AUTO_TEST_CASE( test_time_ordered_design_matrix_uses_active_mixed_dimension_vector_data )
{
    auto dataset = std::make_shared< ObservationDataset<> >( );
    const auto link = createOneWayLinkDefinition( "Station1" );
    dataset->addObservationSet(
            angular_position, link, { Eigen::Vector2d( 1.0, 2.0 ), Eigen::Vector2d( 3.0, 4.0 ) }, { 1.0, 2.0 }, receiver );
    dataset->addObservationSet( one_way_range, link, { Eigen::Vector1d::Constant( 5.0 ) }, { 3.0 }, receiver );
    dataset->rejectObservations( ObservationSelectionCondition<>::timeBounds( 2.0, 2.0 ) );
    // Legacy observable ordering puts range first; time ordering puts the angular components first.
    Eigen::Matrix< double, 3, 1 > design;
    design << 50.0, 10.0, 20.0;
    std::vector< int > timeOrder;
    const auto sorted = simulation_setup::getTimeOrderedDesignMatrix( dataset, design, timeOrder );
    Eigen::Vector3d expected( 10.0, 20.0, 50.0 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( sorted.first, expected, 1.0E-15 );
    const std::vector< double > expectedTimes = { 1.0, 1.0, 3.0 };
    BOOST_CHECK_EQUAL_COLLECTIONS( sorted.second.begin( ), sorted.second.end( ), expectedTimes.begin( ), expectedTimes.end( ) );
}

//! Verify that the covariance as a function of time uses the full weights matrix, with rows and columns ordered by time.
BOOST_AUTO_TEST_CASE( test_covariance_history_with_off_diagonal_weights )
{
    auto dataset = std::make_shared< ObservationDataset<> >( );
    dataset->addObservationSet( one_way_range,
                                createOneWayLinkDefinition( "Station1" ),
                                { Eigen::Vector1d::Zero( ), Eigen::Vector1d::Zero( ), Eigen::Vector1d::Zero( ) },
                                { 3.0, 1.0, 2.0 },
                                receiver );
    Eigen::Matrix3d weights;
    weights << 4.0, 1.0, 0.5, 1.0, 3.0, 0.25, 0.5, 0.25, 2.0;
    const Eigen::SparseMatrix< double > sparseWeights = weights.sparseView( );
    const Eigen::VectorXd diagonal = weights.diagonal( );
    const Eigen::MatrixXd design = Eigen::Vector3d( 1.0, 2.0, 3.0 );
    const Eigen::VectorXd normalization = Eigen::Vector1d::Constant( 2.0 );
    const Eigen::MatrixXd prior = Eigen::MatrixXd::Identity( 1, 1 );

    // Expected covariance is 1 / ( 4 * H^T W H + 1 ), with H and W restricted to the observations up to the output time
    const auto history = simulation_setup::calculateCovarianceUsingDataUpToEpoch(
            dataset, design, normalization, std::vector< double >{ 1.0, 2.0, 3.0 }, sparseWeights, prior );
    BOOST_CHECK_CLOSE_FRACTION( history.at( 1.0 )( 0, 0 ), 1.0 / 49.0, 1.0E-14 );
    BOOST_CHECK_CLOSE_FRACTION( history.at( 2.0 )( 0, 0 ), 1.0 / 133.0, 1.0E-14 );
    BOOST_CHECK_CLOSE_FRACTION( history.at( 3.0 )( 0, 0 ), 1.0 / 177.0, 1.0E-14 );

    auto input = std::make_shared< simulation_setup::EstimationInput<> >( dataset, prior );
    auto output = std::make_shared< simulation_setup::EstimationOutput<> >( Eigen::Vector1d::Zero( ),
                                                                            Eigen::Vector3d::Zero( ),
                                                                            design,
                                                                            diagonal,
                                                                            normalization,
                                                                            Eigen::MatrixXd::Constant( 1, 1, 44.25 ),
                                                                            0.0,
                                                                            0 );
    output->weightsMatrix_ = sparseWeights;
    const auto outputHistory = simulation_setup::calculateCovarianceUsingDataUpToEpoch( input, output, std::vector< double >{ 3.0 } );
    const auto timeStepHistory = simulation_setup::calculateCovarianceUsingDataUpToEpoch( input, output, 1.0 );
    BOOST_CHECK_CLOSE_FRACTION( outputHistory.at( 3.0 )( 0, 0 ), 1.0 / 177.0, 1.0E-14 );
    BOOST_CHECK_CLOSE_FRACTION( timeStepHistory.at( 2.0 )( 0, 0 ), 1.0 / 133.0, 1.0E-14 );
}

//! Verify inspection getters preserve ordering, alignment, and weight correlations.
BOOST_AUTO_TEST_CASE( test_inspection_ordering_alignment_and_correlations )
{
    ObservationDataset<> dataset;
    const auto link1 = createOneWayLinkDefinition( "Station1" );
    const auto link2 = createOneWayLinkDefinition( "Station2" );
    dataset.addObservationSet( angular_position, link2, { Eigen::Vector2d( 10, 11 ), Eigen::Vector2d( 20, 21 ) }, { 4.0, 4.0 }, receiver );
    dataset.addObservationSet(
            one_way_range, link2, { Eigen::Vector1d::Constant( 30 ), Eigen::Vector1d::Constant( 40 ) }, { 3.0, 1.0 }, receiver );
    dataset.addObservationSet( one_way_range, link1, { Eigen::Vector1d::Constant( 50 ) }, { 4.0 }, receiver );
    dataset.addObservationsToSet( 0, { Eigen::Vector2d( 60, 61 ) }, { 2.0 }, {}, {}, {}, true );
    for( unsigned int set = 0; set < 3; ++set )
    {
        auto values = dataset.getObservationsForSet( set );
        for( auto& value : values )
        {
            value *= 0.01;
        }
        dataset.setResidualsForSet( set, values );
        std::vector< Eigen::VectorXd > dependent;
        for( const auto id : dataset.getObservationIdsForSet( set ) )
        {
            dependent.push_back( Eigen::Vector1d::Constant( 100 + id ) );
        }
        dataset.setDependentVariablesForSet( set, dependent );
    }
    Eigen::MatrixXd weights = Eigen::MatrixXd::Identity( 9, 9 );
    for( int i = 0; i < 9; ++i )
    {
        weights( i, i ) = 10 + i;
        for( int j = i + 1; j < 9; ++j )
        {
            weights( i, j ) = weights( j, i ) = 0.001 * ( 10 * i + j );
        }
    }
    dataset.setWeightBlock( { 0, 1, 2, 3, 4, 5 }, { 0, 1, 2, 3, 4, 5 }, weights );
    dataset.rejectObservations( ObservationSelectionCondition<>::timeBounds( 4.0, 4.0 ), "equal epochs" );
    const auto originalRows = dataset.getRows( );
    const auto originalObservations = dataset.getObservations( );
    const auto all = ObservationSelectionCondition<>::all( );
    const ObservationSelectionCondition<> noncontiguous( []( const auto&, const int id ) { return id != 1 && id != 3; } );
    checkIds( dataset.getObservationIds( ), { 0, 1, 2, 3, 4, 5 } );
    checkIds( dataset.getObservationIds( all, ObservationOrdering::estimation ), { 4, 2, 3, 5, 0, 1 } );
    for( const auto& condition : { all, noncontiguous, ObservationSelectionCondition<>::active( ) } )
    {
        for( const auto ordering : { ObservationOrdering::internal, ObservationOrdering::estimation } )
        {
            const auto ids = dataset.getObservationIds( condition, ordering );
            const auto times = dataset.getTimes( condition, ordering );
            const auto values = dataset.getObservations( condition, ordering );
            const auto residuals = dataset.getResiduals( condition, ordering );
            const auto sets = dataset.getSetIds( condition, ordering );
            const auto rows = dataset.getRows( condition, ordering );
            const auto dependent = dataset.getDependentVariableValues( condition, ordering );
            const auto components = dataset.getScalarComponents( condition, ordering );
            const Eigen::MatrixXd selectedWeights = dataset.getWeightMatrix( condition, ordering );
            const auto diagonal = dataset.getWeightDiagonal( condition, ordering );
            BOOST_REQUIRE_EQUAL( times.size( ), ids.size( ) );
            BOOST_REQUIRE_EQUAL( components.size( ), selectedWeights.rows( ) );
            unsigned int scalar = 0;
            for( std::size_t event = 0; event < ids.size( ); ++event )
            {
                const auto id = ids.at( event );
                BOOST_CHECK_EQUAL( times.at( event ), originalRows.at( id ).time_ );
                BOOST_CHECK_EQUAL( sets.at( event ), originalRows.at( id ).setId_ );
                BOOST_CHECK( values.at( event ) == originalObservations.at( id ) );
                BOOST_CHECK( residuals.at( event ) == 0.01 * values.at( event ) );
                BOOST_CHECK_EQUAL( dependent.at( event )( 0 ), 100 + id );
                for( unsigned int component = 0; component < rows.at( event ).scalarSize_; ++component, ++scalar )
                {
                    BOOST_CHECK_EQUAL( components.at( scalar ).first, id );
                    BOOST_CHECK_EQUAL( components.at( scalar ).second, component );
                }
            }
            for( unsigned int i = 0; i < components.size( ); ++i )
            {
                const auto sourceI = originalRows.at( components.at( i ).first ).firstScalarComponent_ + components.at( i ).second;
                BOOST_CHECK_EQUAL( diagonal( i ), weights( sourceI, sourceI ) );
                for( unsigned int j = 0; j < components.size( ); ++j )
                {
                    const auto sourceJ = originalRows.at( components.at( j ).first ).firstScalarComponent_ + components.at( j ).second;
                    BOOST_CHECK_EQUAL( selectedWeights( i, j ), weights( sourceI, sourceJ ) );
                }
            }
            if( ordering == ObservationOrdering::estimation )
            {
                const auto vectorData = dataset.createNewAndKeep( condition )->createObservationVectorData( true );
                BOOST_CHECK( concatenateInspectionValues( values ) == vectorData.getObservationVector( ) );
                BOOST_CHECK( concatenateInspectionValues( residuals ) == vectorData.getResidualVector( ) );
                BOOST_CHECK( selectedWeights == vectorData.getSparseWeightMatrix( ).toDense( ) );
                for( unsigned int i = 0; i < components.size( ); ++i )
                {
                    BOOST_CHECK_EQUAL( vectorData.getObservationIds( ).at( i ), components.at( i ).first );
                    BOOST_CHECK_EQUAL( vectorData.getTimes( ).at( i ), dataset.getObservationTime( components.at( i ).first ) );
                }
            }
        }
        auto internal = dataset.getObservationIds( condition );
        auto estimation = dataset.getObservationIds( condition, ObservationOrdering::estimation );
        std::sort( estimation.begin( ), estimation.end( ) );
        BOOST_CHECK( internal == estimation );
    }
    BOOST_CHECK( dataset.getObservations( ) == originalObservations );
}

//! Verify detached inspection snapshots survive dataset mutation and destruction.
BOOST_AUTO_TEST_CASE( test_inspection_snapshots_survive_mutation_and_destruction )
{
    ObservationDataset<>::InspectionMetadata metadata;
    std::vector< Eigen::VectorXd > values, residuals, dependent;
    std::vector< ObservationDatasetRow<> > rows;
    Eigen::VectorXd diagonal;
    Eigen::SparseMatrix< double > weights;
    {
        ObservationDataset<> dataset;
        const auto link = createOneWayLinkDefinition( "Station1" );
        auto bookkeeping = std::make_shared< simulation_setup::ObservationDependentVariableBookkeeping >( angular_position, link );
        bookkeeping->addDependentVariable( simulation_setup::elevationAngleDependentVariable( receiver, LinkEndId( "Vehicle", "" ) ) );
        auto ancillary = std::make_shared< ObservationAncillarySimulationSettings >( );
        ancillary->setAncillaryDoubleData( doppler_integration_time, 10 );
        dataset.addObservationSet( angular_position,
                                   link,
                                   { Eigen::Vector2d( 1, 2 ), Eigen::Vector2d( 3, 4 ) },
                                   { 2.0, 1.0 },
                                   receiver,
                                   { Eigen::Vector1d::Constant( 5 ), Eigen::Vector1d::Constant( 6 ) },
                                   bookkeeping,
                                   ancillary,
                                   {},
                                   { Eigen::Vector2d( 0.1, 0.2 ), Eigen::Vector2d( 0.3, 0.4 ) } );
        Eigen::Matrix4d initialWeights = Eigen::Matrix4d::Identity( );
        initialWeights( 0, 3 ) = initialWeights( 3, 0 ) = 0.25;
        dataset.setWeightMatrixForSet( 0, initialWeights );
        values = dataset.getObservations( );
        residuals = dataset.getResiduals( );
        dependent = dataset.getDependentVariableValues( );
        rows = dataset.getRows( );
        weights = dataset.getWeightMatrix( );
        diagonal = dataset.getWeightDiagonal( );
        metadata = dataset.getMetadata( );
        dataset.setObservationsForSet( 0, { Eigen::Vector2d( 11, 12 ), Eigen::Vector2d( 13, 14 ) } );
        dataset.setResidualsForSet( 0, { Eigen::Vector2d( 1.1, 1.2 ), Eigen::Vector2d( 1.3, 1.4 ) } );
        dataset.setDependentVariablesForSet( 0, { Eigen::Vector1d::Constant( 15 ), Eigen::Vector1d::Constant( 16 ) } );
        dataset.setConstantSingleObservationScalarWeight( ObservationSelectionCondition<>::setId( 0 ), 4 );
        dataset.setWeightBlock( { 0 }, { 1 }, Eigen::Matrix2d::Constant( 0.5 ) );
        dataset.getAncillarySettingsForSet( 0 )->setAncillaryDoubleData( doppler_integration_time, 20 );
        dataset.getDependentVariableBookkeeping( dataset.getObservationSetMetadata( 0 ).dependentVariableLayoutId_ )
                ->getDependentVariableSettings( )
                .front( )
                ->linkEndId_ = LinkEndId( "Moon", "" );
        dataset.setLinkEndReferencePoint( "Earth", "OtherStation", transmitter );
        dataset.rejectObservations( ObservationSelectionCondition<>::all( ), "test" );
        dataset.restoreObservations( ObservationSelectionCondition<>::all( ) );
        dataset.addObservationsToSet( 0, {}, {}, {}, {}, {}, true );
        dataset.addObservationsToSet( 0, { Eigen::Vector2d( 7, 8 ) }, { 0.0 }, { Eigen::Vector1d::Constant( 17 ) } );
        dataset.removeObservations( ObservationSelectionCondition<>::timeBounds( 1.0, 1.0 ) );
        BOOST_CHECK_EQUAL( dataset.getObservations( ).front( )( 0 ), 11 );
        BOOST_CHECK_EQUAL( values.front( )( 0 ), 1 );
        BOOST_CHECK_EQUAL( residuals.front( )( 0 ), 0.1 );
        BOOST_CHECK_EQUAL( dependent.front( )( 0 ), 5 );
        BOOST_CHECK_EQUAL( weights.coeff( 0, 3 ), 0.25 );
        BOOST_CHECK_EQUAL( diagonal( 0 ), 1 );
        BOOST_CHECK_EQUAL( rows.at( 1 ).observationId_, 1 );
        BOOST_CHECK_EQUAL( rows.at( 0 ).time_, 2 );
        BOOST_CHECK( rows.at( 0 ).isActive_ );
        BOOST_CHECK( rows.at( 0 ).rejectionReason_.empty( ) );
        BOOST_CHECK_EQUAL( std::get< 2 >( metadata.at( 0 ) )->getAncillaryDoubleData( doppler_integration_time ), 10 );
        BOOST_CHECK_EQUAL( std::get< 3 >( metadata.at( 0 ) ).begin( )->second->linkEndId_.bodyName_, "Vehicle" );
        BOOST_CHECK( std::get< 1 >( metadata.at( 0 ) ) == link );
        // Change returned nested objects; every corresponding source field must remain independent.
        values.front( )( 0 ) = -1;
        residuals.front( )( 0 ) = -2;
        dependent.front( )( 0 ) = -3;
        rows.front( ).dependentVariableValues_( 0 ) = -4;
        diagonal( 0 ) = -5;
        weights.coeffRef( 0, 0 ) = -6;
        std::get< 2 >( metadata.at( 0 ) )->setAncillaryDoubleData( doppler_integration_time, 30 );
        std::get< 3 >( metadata.at( 0 ) ).begin( )->second->linkEndId_ = LinkEndId( "Mars", "" );
        BOOST_CHECK_EQUAL( dataset.getObservationValue( 0 )( 0 ), 11 );
        BOOST_CHECK_EQUAL( dataset.getResidualValue( 0 )( 0 ), 1.1 );
        BOOST_CHECK_EQUAL( dataset.getDependentVariables( 0 )( 0 ), 15 );
        BOOST_CHECK_EQUAL( dataset.getWeightValue( 0 )( 0 ), 4 );
        BOOST_CHECK_EQUAL( dataset.getAncillarySettingsForSet( 0 )->getAncillaryDoubleData( doppler_integration_time ), 20 );
        BOOST_CHECK_EQUAL( dataset.getDependentVariableBookkeeping( dataset.getObservationSetMetadata( 0 ).dependentVariableLayoutId_ )
                                   ->getDependentVariableSettings( )
                                   .front( )
                                   ->linkEndId_.bodyName_,
                           "Moon" );
    }
    BOOST_CHECK_EQUAL( values.at( 1 )( 1 ), 4 );
    BOOST_CHECK_EQUAL( residuals.at( 1 )( 1 ), 0.4 );
    BOOST_CHECK_EQUAL( dependent.at( 1 )( 0 ), 6 );
    BOOST_CHECK_EQUAL( rows.at( 1 ).time_, 1 );
    BOOST_CHECK_EQUAL( weights.coeff( 3, 0 ), 0.25 );
    BOOST_CHECK_EQUAL( std::get< 2 >( metadata.at( 0 ) )->getAncillaryDoubleData( doppler_integration_time ), 30 );
}

//! Verify inspection handles empty data, missing fields, and precise times.
BOOST_AUTO_TEST_CASE( test_inspection_empty_missing_and_precise_times )
{
    ObservationDataset<> empty;
    BOOST_CHECK( empty.getTimes( ).empty( ) );
    BOOST_CHECK( empty.getObservations( ).empty( ) );
    BOOST_CHECK( empty.getResiduals( ).empty( ) );
    BOOST_CHECK( empty.getRows( ).empty( ) );
    BOOST_CHECK( empty.getMetadata( ).empty( ) );
    BOOST_CHECK_EQUAL( empty.getWeightMatrix( ).rows( ), 0 );
    BOOST_CHECK_EQUAL( empty.getWeightDiagonal( ).size( ), 0 );
    const auto link = createOneWayLinkDefinition( "Station1" );
    empty.addObservationSet( one_way_range, link, { Eigen::Vector1d::Constant( 1 ) }, { 2.0 }, receiver );
    BOOST_CHECK_EQUAL( empty.getResiduals( ).front( )( 0 ), 0 );  // Existing storage convention for omitted residuals.
    BOOST_CHECK_EQUAL( empty.getDependentVariableValues( ).front( ).size( ), 0 );
    BOOST_CHECK( empty.getTimes( !ObservationSelectionCondition<>::all( ) ).empty( ) );

    using Scalar = long double;
    using Vector = Eigen::Matrix< Scalar, Eigen::Dynamic, 1 >;
    ObservationDataset< Scalar, Time > precise;
    const Time t0( 1000000, 0.000000001L ), t1( 1000000, 0.000000002L );
    const Scalar x = std::nextafter( 1.0L, 2.0L );
    precise.addObservationSet( one_way_range,
                               link,
                               { Vector::Constant( 1, x ), Vector::Constant( 1, 2 * x ) },
                               { t1, t0 },
                               receiver,
                               {},
                               nullptr,
                               nullptr,
                               {},
                               { Vector::Constant( 1, x ), Vector::Constant( 1, -x ) } );
    BOOST_CHECK_EQUAL( static_cast< double >( t0 ), static_cast< double >( t1 ) );
    BOOST_CHECK( t0 != t1 );
    for( const auto order : { ObservationOrdering::internal, ObservationOrdering::estimation } )
    {
        const auto times = precise.getTimes( ObservationSelectionCondition< Scalar, Time >::all( ), order );
        BOOST_CHECK( times.at( 0 ) == t1 );
        BOOST_CHECK( times.at( 1 ) == t0 );
        BOOST_CHECK( precise.getObservations( ObservationSelectionCondition< Scalar, Time >::all( ), order ).front( )( 0 ) == x );
        BOOST_CHECK( precise.getResiduals( ObservationSelectionCondition< Scalar, Time >::all( ), order ).front( )( 0 ) == x );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
