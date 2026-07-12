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
#include <iostream>
#include <typeinfo>

#include <boost/test/unit_test.hpp>
#include <Eigen/Core>

#include <cereal/archives/binary.hpp>

#include "tudat/io/serialization.h"

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
                                             const std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > >& rawSolution ) {
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

std::shared_ptr< propagators::SingleArcVariationalSimulationResults< double, double > > createSingleArcVariationalSimulationResults( )
{
    auto singleArcResults = createSingleArcSimulationResults( );
    auto variationalResults =
            std::make_shared< propagators::SingleArcVariationalSimulationResults< double, double > >( singleArcResults, 2, 1 );

    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > > fullSolution;
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > solutionAtT0( 2, 4 );
    solutionAtT0 << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0;
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > solutionAtT1( 2, 4 );
    solutionAtT1 << 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0;
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

std::shared_ptr< propagators::MultiArcSimulationResults< propagators::SingleArcSimulationResults, double, double > >
createMultiArcSimulationResults( )
{
    std::vector< std::shared_ptr< propagators::SingleArcSimulationResults< double, double > > > singleArcResults = {
        createSingleArcSimulationResults( ), createSingleArcSimulationResults( )
    };

    std::vector< std::shared_ptr< propagators::SingleArcDependentVariablesInterface< double > > > singleArcInterfaces = {
        singleArcResults.at( 0 )->getSingleArcDependentVariablesInterface( ),
        singleArcResults.at( 1 )->getSingleArcDependentVariablesInterface( )
    };

    auto dependentVariableInterface = std::make_shared< propagators::MultiArcDependentVariablesInterface< double > >(
            singleArcInterfaces, std::vector< double >( ), std::vector< double >( ) );

    auto results = std::make_shared< propagators::MultiArcSimulationResults< propagators::SingleArcSimulationResults, double, double > >(
            singleArcResults, dependentVariableInterface );
    results->setPropagationIsPerformed( );
    return results;
}

std::shared_ptr< propagators::HybridArcSimulationResults< propagators::SingleArcSimulationResults, double, double > >
createHybridArcSimulationResults( )
{
    return std::make_shared< propagators::HybridArcSimulationResults< propagators::SingleArcSimulationResults, double, double > >(
            createSingleArcSimulationResults( ), createMultiArcSimulationResults( ) );
}

template< typename ResultType >
std::shared_ptr< ResultType > roundTripSerialize( const std::shared_ptr< ResultType >& originalResult )
{
    std::stringstream serializationStream;

    {
        cereal::BinaryOutputArchive outputArchive( serializationStream );
        outputArchive( originalResult );
    }

    std::shared_ptr< ResultType > deserializedResult;

    {
        cereal::BinaryInputArchive inputArchive( serializationStream );
        inputArchive( deserializedResult );
    }

    return deserializedResult;
}

}  // namespace

BOOST_AUTO_TEST_SUITE( test_PropagationResults_serialization )

BOOST_AUTO_TEST_CASE( test_PropagationResultsSerialization )
{
    std::vector< std::shared_ptr< propagators::SimulationResults< double, double > > > resultsVector = {
        createSingleArcSimulationResults( ),
        createSingleArcVariationalSimulationResults( ),
        createMultiArcSimulationResults( ),
        createHybridArcSimulationResults( )
    };

    for( const auto& results : resultsVector )
    {
        auto deserializedResults = roundTripSerialize( results );

        BOOST_REQUIRE( deserializedResults != nullptr );
        if( !( *results == *deserializedResults ) )
        {
            std::cout << "--- Serialization mismatch for type: " << typeid( *results ).name( ) << " ---\n";
            // Single-arc
            if( auto single = std::dynamic_pointer_cast< propagators::SingleArcSimulationResults< double, double > >( results ) )
            {
                auto sdes = std::dynamic_pointer_cast< propagators::SingleArcSimulationResults< double, double > >( deserializedResults );
                if( sdes )
                {
                    std::cout << "single: rawSolution size (orig/deser): " << single->getEquationsOfMotionNumericalSolutionRaw( ).size( )
                              << " / " << sdes->getEquationsOfMotionNumericalSolutionRaw( ).size( ) << "\n";
                    std::cout << "single: dependentVar size (orig/deser): " << single->getDependentVariableHistoryDouble( ).size( ) << " / "
                              << sdes->getDependentVariableHistoryDouble( ).size( ) << "\n";
                    std::cout << "single: compTime size (orig/deser): " << single->getCumulativeComputationTimeHistory( ).size( ) << " / "
                              << sdes->getCumulativeComputationTimeHistory( ).size( ) << "\n";
                    std::cout << "single: funcEval size (orig/deser): " << single->getCumulativeNumberOfFunctionEvaluations( ).size( )
                              << " / " << sdes->getCumulativeNumberOfFunctionEvaluations( ).size( ) << "\n";
                }
            }
            // Variational
            if( auto var = std::dynamic_pointer_cast< propagators::SingleArcVariationalSimulationResults< double, double > >( results ) )
            {
                auto vdes = std::dynamic_pointer_cast< propagators::SingleArcVariationalSimulationResults< double, double > >(
                        deserializedResults );
                if( vdes )
                {
                    std::cout << "variational: stateTrans size (orig/deser): " << var->getStateTransitionSolution( ).size( ) << " / "
                              << vdes->getStateTransitionSolution( ).size( ) << "\n";
                    std::cout << "variational: sens size (orig/deser): " << var->getSensitivitySolution( ).size( ) << " / "
                              << vdes->getSensitivitySolution( ).size( ) << "\n";
                    auto dyn = var->getDynamicsResults( );
                    auto ddyn = vdes->getDynamicsResults( );
                    if( dyn && ddyn )
                    {
                        std::cout << "variational->dynamics: rawSolution size (orig/deser): "
                                  << dyn->getEquationsOfMotionNumericalSolutionRaw( ).size( ) << " / "
                                  << ddyn->getEquationsOfMotionNumericalSolutionRaw( ).size( ) << "\n";
                    }
                }
            }
            // Multi-arc
            if( auto multi = std::dynamic_pointer_cast<
                        propagators::MultiArcSimulationResults< propagators::SingleArcSimulationResults, double, double > >( results ) )
            {
                auto mdes = std::dynamic_pointer_cast<
                        propagators::MultiArcSimulationResults< propagators::SingleArcSimulationResults, double, double > >(
                        deserializedResults );
                if( mdes )
                {
                    std::cout << "multi: single-arc count (orig/deser): " << multi->getSingleArcResults( ).size( ) << " / "
                              << mdes->getSingleArcResults( ).size( ) << "\n";
                }
            }
            // Hybrid
            if( auto hyb = std::dynamic_pointer_cast<
                        propagators::HybridArcSimulationResults< propagators::SingleArcSimulationResults, double, double > >( results ) )
            {
                auto hdes = std::dynamic_pointer_cast<
                        propagators::HybridArcSimulationResults< propagators::SingleArcSimulationResults, double, double > >(
                        deserializedResults );
                if( hdes )
                {
                    auto s = hyb->getSingleArcResults( );
                    auto sd = hdes->getSingleArcResults( );
                    std::cout << "hybrid: single raw size (orig/deser): " << s->getEquationsOfMotionNumericalSolutionRaw( ).size( ) << " / "
                              << sd->getEquationsOfMotionNumericalSolutionRaw( ).size( ) << "\n";
                }
            }
        }
        BOOST_CHECK( *results == *deserializedResults );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat