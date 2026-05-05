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

#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <tuple>
#include <vector>

#include <boost/test/unit_test.hpp>
#include <Eigen/Core>

#include <cereal/archives/binary.hpp>

#include "tudat/io/serialization/base.h"

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/propagation_setup/dependentVariablesInterface.h"
#include "tudat/simulation/propagation_setup/propagationProcessingSettings.h"
#include "tudat/simulation/propagation_setup/propagationResults.h"
#include "tudat/simulation/propagation_setup/propagationTermination.h"

namespace tudat
{
namespace unit_tests
{

namespace
{

template< typename MapType >
void compareMaps( const MapType& left, const MapType& right )
{
    BOOST_REQUIRE_EQUAL( left.size( ), right.size( ) );
    auto leftIterator = left.begin( );
    auto rightIterator = right.begin( );
    for( ; leftIterator != left.end( ); ++leftIterator, ++rightIterator )
    {
        BOOST_CHECK_EQUAL( leftIterator->first, rightIterator->first );
    }
}

std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > createRawStateHistory( )
{
    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > rawStateHistory;

    Eigen::Vector2d stateAtT0;
    stateAtT0 << 1.0, 2.0;
    rawStateHistory[ 0.0 ] = stateAtT0;

    Eigen::Vector2d stateAtT1;
    stateAtT1 << 3.0, 4.0;
    rawStateHistory[ 10.0 ] = stateAtT1;

    return rawStateHistory;
}

std::map< double, Eigen::VectorXd > createDependentVariableHistory( )
{
    std::map< double, Eigen::VectorXd > dependentVariableHistory;
    dependentVariableHistory[ 0.0 ] = ( Eigen::Vector2d( ) << 5.0, 6.0 ).finished( );
    dependentVariableHistory[ 10.0 ] = ( Eigen::Vector2d( ) << 7.0, 8.0 ).finished( );
    return dependentVariableHistory;
}

std::map< double, double > createComputationTimeHistory( )
{
    return { { 0.0, 1.5 }, { 10.0, 2.5 } };
}

std::map< double, unsigned int > createFunctionEvaluationHistory( )
{
    return { { 0.0, 12 }, { 10.0, 24 } };
}

std::shared_ptr< propagators::SingleArcSimulationResults< double, double > > createSingleArcSimulationResults( )
{
        std::map< tudat::propagators::IntegratedStateType,
                            std::vector< std::tuple< std::string, std::string, tudat::propagators::PropagatorType > > >
                        integratedStateAndBodyList;

    auto rawSolutionConversionFunction = []( std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > >& processedSolution,
                                             const std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > >& rawSolution )
    {
        processedSolution = rawSolution;
    };

    auto results = std::make_shared< propagators::SingleArcSimulationResults< double, double > >(
            integratedStateAndBodyList,
            std::shared_ptr< propagators::SingleArcPropagatorProcessingSettings >( ),
            rawSolutionConversionFunction,
            std::shared_ptr< propagators::SingleArcDependentVariablesInterface< double > >( ),
            true );

    results->reset( createRawStateHistory( ),
                    createDependentVariableHistory( ),
                    createComputationTimeHistory( ),
                    createFunctionEvaluationHistory( ),
                    std::make_shared< propagators::PropagationTerminationDetails >( propagators::propagation_never_run ) );
    results->finalizePropagation( createFunctionEvaluationHistory( ) );
    return results;
}

void compareSingleArcResults( propagators::SingleArcSimulationResults< double, double >& left,
                              propagators::SingleArcSimulationResults< double, double >& right )
{
    compareMaps( left.getEquationsOfMotionNumericalSolutionRaw( ), right.getEquationsOfMotionNumericalSolutionRaw( ) );
    compareMaps( left.getEquationsOfMotionNumericalSolution( ), right.getEquationsOfMotionNumericalSolution( ) );
    compareMaps( left.getDependentVariableHistory( ), right.getDependentVariableHistory( ) );
    compareMaps( left.getCumulativeComputationTimeHistoryTimeType( ), right.getCumulativeComputationTimeHistoryTimeType( ) );
    compareMaps( left.getCumulativeNumberOfFunctionEvaluationsTimeType( ), right.getCumulativeNumberOfFunctionEvaluationsTimeType( ) );
    BOOST_CHECK_EQUAL( left.getPropagationIsPerformed( ), right.getPropagationIsPerformed( ) );
    BOOST_CHECK_EQUAL( left.getSolutionIsCleared( ), right.getSolutionIsCleared( ) );
    BOOST_CHECK_EQUAL( left.getPropagationTerminationReason( )->getPropagationTerminationReason( ),
                       right.getPropagationTerminationReason( )->getPropagationTerminationReason( ) );
    BOOST_CHECK_EQUAL( left.getPropagationTerminationReason( )->getTerminationOnExactCondition( ),
                       right.getPropagationTerminationReason( )->getTerminationOnExactCondition( ) );
}

std::shared_ptr< propagators::SingleArcVariationalSimulationResults< double, double > > createSingleArcVariationalResults( )
{
    std::shared_ptr< propagators::SingleArcSimulationResults< double, double > > singleArcResults = createSingleArcSimulationResults( );
    auto variationalResults = std::make_shared< propagators::SingleArcVariationalSimulationResults< double, double > >(
            singleArcResults, 2, 1 );

    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > > fullSolution;
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > solutionAtT0( 2, 4 );
    solutionAtT0 << 1.0, 2.0, 3.0, 4.0,
                    5.0, 6.0, 7.0, 8.0;
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > solutionAtT1( 2, 4 );
    solutionAtT1 << 9.0, 10.0, 11.0, 12.0,
                    13.0, 14.0, 15.0, 16.0;
    fullSolution[ 0.0 ] = solutionAtT0;
    fullSolution[ 10.0 ] = solutionAtT1;

    variationalResults->reset( fullSolution,
                               createDependentVariableHistory( ),
                               createComputationTimeHistory( ),
                               createFunctionEvaluationHistory( ),
                               std::make_shared< propagators::PropagationTerminationDetails >( propagators::propagation_never_run ) );
    variationalResults->finalizePropagation( createFunctionEvaluationHistory( ) );
    return variationalResults;
}

void compareVariationalResults( propagators::SingleArcVariationalSimulationResults< double, double >& left,
                                propagators::SingleArcVariationalSimulationResults< double, double >& right )
{
    compareSingleArcResults( *left.getDynamicsResults( ), *right.getDynamicsResults( ) );
    compareMaps( left.getStateTransitionSolution( ), right.getStateTransitionSolution( ) );
    compareMaps( left.getSensitivitySolution( ), right.getSensitivitySolution( ) );
    BOOST_CHECK_EQUAL( left.getStateTransitionMatrixSize( ), right.getStateTransitionMatrixSize( ) );
    BOOST_CHECK_EQUAL( left.getSensitivityMatrixSize( ), right.getSensitivityMatrixSize( ) );
}

std::shared_ptr< propagators::MultiArcSimulationResults< propagators::SingleArcSimulationResults, double, double > >
createMultiArcResults( )
{
    std::vector< std::shared_ptr< propagators::SingleArcDependentVariablesInterface< double > > > singleArcInterfaces;
    std::vector< double > arcStartTimes;
    std::vector< double > arcEndTimes;

    auto dependentVariableInterface = std::make_shared< propagators::MultiArcDependentVariablesInterface< double > >(
            singleArcInterfaces, arcStartTimes, arcEndTimes );

    auto results = std::make_shared< propagators::MultiArcSimulationResults< propagators::SingleArcSimulationResults, double, double > >(
            std::vector< std::shared_ptr< propagators::SingleArcSimulationResults< double, double > > >( ),
            dependentVariableInterface );
    results->setPropagationIsPerformed( );
    return results;
}

std::shared_ptr< propagators::HybridArcSimulationResults< propagators::SingleArcSimulationResults, double, double > >
createHybridArcResults( )
{
    auto singleArcResults = createSingleArcSimulationResults( );
    auto multiArcResults = createMultiArcResults( );
    return std::make_shared< propagators::HybridArcSimulationResults< propagators::SingleArcSimulationResults, double, double > >(
            singleArcResults, multiArcResults );
}

}  // namespace

BOOST_AUTO_TEST_SUITE( test_PropagationResults_serialization )

BOOST_AUTO_TEST_CASE( test_SingleArcSimulationResultsSerialization )
{
    auto originalResults = createSingleArcSimulationResults( );

    std::stringstream serializationStream;
    {
        cereal::BinaryOutputArchive outputArchive( serializationStream );
        outputArchive( originalResults );
    }

    std::shared_ptr< propagators::SingleArcSimulationResults< double, double > > deserializedResults;
    {
        cereal::BinaryInputArchive inputArchive( serializationStream );
        inputArchive( deserializedResults );
    }

    BOOST_REQUIRE( deserializedResults != nullptr );
    compareSingleArcResults( *originalResults, *deserializedResults );
}

BOOST_AUTO_TEST_CASE( test_SingleArcVariationalSimulationResultsSerialization )
{
    auto originalResults = createSingleArcVariationalResults( );

    std::stringstream serializationStream;
    {
        cereal::BinaryOutputArchive outputArchive( serializationStream );
        outputArchive( originalResults );
    }

    std::shared_ptr< propagators::SingleArcVariationalSimulationResults< double, double > > deserializedResults;
    {
        cereal::BinaryInputArchive inputArchive( serializationStream );
        inputArchive( deserializedResults );
    }

    BOOST_REQUIRE( deserializedResults != nullptr );
    compareVariationalResults( *originalResults, *deserializedResults );
}

BOOST_AUTO_TEST_CASE( test_MultiArcSimulationResultsSerialization )
{
    auto originalResults = createMultiArcResults( );

    std::stringstream serializationStream;
    {
        cereal::BinaryOutputArchive outputArchive( serializationStream );
        outputArchive( originalResults );
    }

    std::shared_ptr< propagators::MultiArcSimulationResults< propagators::SingleArcSimulationResults, double, double > > deserializedResults;
    {
        cereal::BinaryInputArchive inputArchive( serializationStream );
        inputArchive( deserializedResults );
    }

    BOOST_REQUIRE( deserializedResults != nullptr );
    BOOST_CHECK_EQUAL( originalResults->getPropagationIsPerformed( ), deserializedResults->getPropagationIsPerformed( ) );
    BOOST_CHECK_EQUAL( originalResults->getSolutionIsCleared( ), deserializedResults->getSolutionIsCleared( ) );
    BOOST_CHECK_EQUAL( originalResults->getSingleArcResults( ).size( ), deserializedResults->getSingleArcResults( ).size( ) );
    BOOST_CHECK_EQUAL( originalResults->getArcStartTimes( ).size( ), deserializedResults->getArcStartTimes( ).size( ) );
    BOOST_CHECK_EQUAL( originalResults->getArcEndTimes( ).size( ), deserializedResults->getArcEndTimes( ).size( ) );
}

BOOST_AUTO_TEST_CASE( test_HybridArcSimulationResultsSerialization )
{
    auto originalResults = createHybridArcResults( );

    std::stringstream serializationStream;
    {
        cereal::BinaryOutputArchive outputArchive( serializationStream );
        outputArchive( originalResults );
    }

    std::shared_ptr< propagators::HybridArcSimulationResults< propagators::SingleArcSimulationResults, double, double > > deserializedResults;
    {
        cereal::BinaryInputArchive inputArchive( serializationStream );
        inputArchive( deserializedResults );
    }

    BOOST_REQUIRE( deserializedResults != nullptr );
    BOOST_CHECK_EQUAL( originalResults->getSingleArcResults( )->getPropagationIsPerformed( ),
                       deserializedResults->getSingleArcResults( )->getPropagationIsPerformed( ) );
    BOOST_CHECK_EQUAL( originalResults->getMultiArcResults( )->getPropagationIsPerformed( ),
                       deserializedResults->getMultiArcResults( )->getPropagationIsPerformed( ) );
    BOOST_CHECK_EQUAL( originalResults->integrationCompletedSuccessfully( ), deserializedResults->integrationCompletedSuccessfully( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat