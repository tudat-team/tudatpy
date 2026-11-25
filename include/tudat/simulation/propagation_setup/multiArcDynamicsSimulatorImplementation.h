/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_MULTIARCDYNAMICSSIMULATORIMPLEMENTATION_H
#define TUDAT_MULTIARCDYNAMICSSIMULATORIMPLEMENTATION_H

#include "tudat/simulation/propagation_setup/multiArcDynamicsSimulator.h"

namespace tudat
{

namespace propagators
{

template< typename StateScalarType, typename TimeType >
MultiArcDynamicsSimulator< StateScalarType, TimeType >::MultiArcDynamicsSimulator(
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< MultiArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings,
        const bool areEquationsOfMotionToBeIntegrated ):
    DynamicsSimulator< StateScalarType, TimeType >( bodies, propagatorSettings ), multiArcPropagatorSettings_( propagatorSettings )
{
    if( multiArcPropagatorSettings_ == nullptr )
    {
        throw std::runtime_error( "Error when creating multi-arc dynamics simulator, input is not multi arc" );
    }
    else
    {
        std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType, TimeType > > > singleArcSettings =
                multiArcPropagatorSettings_->getSingleArcSettings( );

        // Create dynamics simulators
        std::vector< std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType > > > singleArcResults;
        for( unsigned int i = 0; i < singleArcSettings.size( ); i++ )
        {
            singleArcDynamicsSimulators_.push_back( std::make_shared< SingleArcDynamicsSimulator< StateScalarType, TimeType > >(
                    bodies,
                    singleArcSettings.at( i ),
                    false,
                    PredefinedSingleArcStateDerivativeModels< StateScalarType, TimeType >( ),
                    true ) );
            singleArcResults.push_back( singleArcDynamicsSimulators_.at( i )->getSingleArcPropagationResults( ) );
            if( singleArcDynamicsSimulators_.at( i )->getPropagatorSettings( )->getOutputSettings( )->getSetIntegratedResult( ) )
            {
                singleArcDynamicsSimulators_.at( i )->createAndSetIntegratedStateProcessors( );
            }
        }

        std::vector< std::shared_ptr< SingleArcDependentVariablesInterface< TimeType > > > singleArcInterfaces;
        for( unsigned int i = 0; i < singleArcSettings.size( ); i++ )
        {
            singleArcInterfaces.push_back(
                    singleArcDynamicsSimulators_.at( i )->getSingleArcPropagationResults( )->getSingleArcDependentVariablesInterface( ) );
        }

        std::shared_ptr< MultiArcDependentVariablesInterface< TimeType > > dependentVariableInterface =
                std::make_shared< MultiArcDependentVariablesInterface< TimeType > >(
                        singleArcInterfaces, std::vector< double >( ), std::vector< double >( ) );

        propagationResults_ = std::make_shared< MultiArcResults >( singleArcResults, dependentVariableInterface );

        // Integrate equations of motion if required.
        if( areEquationsOfMotionToBeIntegrated )
        {
            integrateEquationsOfMotion( multiArcPropagatorSettings_->getInitialStates( ) );
        }
    }
}

template< typename StateScalarType, typename TimeType >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > MultiArcDynamicsSimulator< StateScalarType, TimeType >::getArcInitialState(
        const int arcIndex,
        const std::shared_ptr< MultiArcInitialStateProvider< StateScalarType > > stateProvider )
{
    bool initialStateFromPreviousArc = false;
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > currentArcInitialState =
            stateProvider->getArcInitialState( arcIndex, initialStateFromPreviousArc );
    if( initialStateFromPreviousArc )
    {
        currentArcInitialState.block( 0, currentArcInitialState.cols( ) - 1, currentArcInitialState.rows( ), 1 ) =
                getArcInitialStateFromPreviousArcResult(
                        propagationResults_->getSingleArcResults( ).at( arcIndex - 1 )->getEquationsOfMotionNumericalSolution( ),
                        singleArcDynamicsSimulators_.at( arcIndex )->getInitialPropagationTime( ) );
    }
    return currentArcInitialState;
}

template< typename StateScalarType, typename TimeType >
void MultiArcDynamicsSimulator< StateScalarType, TimeType >::integrateEquationsOfMotion(
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& concatenatedInitialStates )
{
    std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > splitInitialState;

    int currentIndex = 0;
    for( unsigned int i = 0; i < singleArcDynamicsSimulators_.size( ); i++ )
    {
        int currentSize = singleArcDynamicsSimulators_.at( i )->getPropagatorSettings( )->getConventionalStateSize( );
        splitInitialState.push_back( concatenatedInitialStates.block( currentIndex, 0, currentSize, 1 ) );
        currentIndex += currentSize;
    }

    if( currentIndex != concatenatedInitialStates.rows( ) )
    {
        throw std::runtime_error( "Error when doing multi-arc integration, input state vector size is incompatible with settings" );
    }

    integrateEquationsOfMotion( splitInitialState );
}

template< typename StateScalarType, typename TimeType >
void MultiArcDynamicsSimulator< StateScalarType, TimeType >::integrate(
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& concatenatedInitialStates )
{
    integrateEquationsOfMotion( concatenatedInitialStates );
}

template< typename StateScalarType, typename TimeType >
void MultiArcDynamicsSimulator< StateScalarType, TimeType >::integrateEquationsOfMotion(
        const std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& initialStatesList )
{
    integrateEquationsOfMotion< MultiArcResults >(
            propagationResults_, std::make_shared< MultiArcInitialStateProvider< StateScalarType > >( initialStatesList ) );
}

template< typename StateScalarType, typename TimeType >
template< typename MultiArcSimulationResults >
void MultiArcDynamicsSimulator< StateScalarType, TimeType >::integrateEquationsOfMotion(
        const std::shared_ptr< MultiArcSimulationResults > propagationResults,
        const std::shared_ptr< MultiArcInitialStateProvider< StateScalarType > > initialStateProvider )
{
    checkPropagationResultsObjectConsistency< StateScalarType, TimeType, MultiArcSimulationResults >( propagationResults_,
                                                                                                      propagationResults );

    Eigen::Matrix< StateScalarType, Eigen::Dynamic, MultiArcSimulationResults::single_arc_type::number_of_columns > currentArcInitialState;
    std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, MultiArcSimulationResults::single_arc_type::number_of_columns > >
            arcInitialStateList;

    initialStateProvider->restartPropagation( );
    propagationResults->restartPropagation( );
    propagationResults_->restartPropagation( );

    printPrePropagationMessages( );

    // Propagate dynamics for each arc
    for( unsigned int i = 0; i < singleArcDynamicsSimulators_.size( ); i++ )
    {
        currentArcInitialState = getArcInitialState( i, initialStateProvider );
        arcInitialStateList.push_back( currentArcInitialState );

        singleArcDynamicsSimulators_.at( i )->template integrateEquationsOfMotion< typename MultiArcSimulationResults::single_arc_type >(
                currentArcInitialState, propagationResults->getSingleArcResults( ).at( i ) );
    }

    printPostPropagationMessages( );
    propagationResults->setPropagationIsPerformed( );
    propagationResults_->setPropagationIsPerformed( );

    if( initialStateProvider->getUpdateInitialStates( ) )
    {
        std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > newInitialStates;
        for( unsigned int i = 0; i < arcInitialStateList.size( ); i++ )
        {
            newInitialStates.push_back( arcInitialStateList.at( i ).block(
                    0, arcInitialStateList.at( i ).cols( ) - 1, arcInitialStateList.at( i ).rows( ), 1 ) );
        }

        multiArcPropagatorSettings_->resetInitialStatesList( newInitialStates );
    }

    processNumericalEquationsOfMotionSolution( );
}

template< typename StateScalarType, typename TimeType >
std::vector< std::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > >
MultiArcDynamicsSimulator< StateScalarType, TimeType >::getDynamicsStateDerivative( )
{
    std::vector< std::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > > dynamicsStateDerivatives;
    for( unsigned int i = 0; i < singleArcDynamicsSimulators_.size( ); i++ )
    {
        dynamicsStateDerivatives.push_back( singleArcDynamicsSimulators_.at( i )->getDynamicsStateDerivative( ) );
    }
    return dynamicsStateDerivatives;
}

template< typename StateScalarType, typename TimeType >
void MultiArcDynamicsSimulator< StateScalarType, TimeType >::processNumericalEquationsOfMotionSolution( )
{
    if( multiArcPropagatorSettings_->getOutputSettings( )->getSetIntegratedResult( ) )
    {
        try
        {
            std::map< IntegratedStateType,
                      std::vector< std::shared_ptr< SingleArcIntegratedStateProcessor< TimeType, StateScalarType > > > >
                    singleArcIntegratedStatesProcessors;

            for( unsigned int i = 0; i < singleArcDynamicsSimulators_.size( ); i++ )
            {
                std::map< IntegratedStateType, std::shared_ptr< SingleArcIntegratedStateProcessor< TimeType, StateScalarType > > >
                        currentArcStateProcessors = singleArcDynamicsSimulators_.at( i )->getIntegratedStateProcessors( );
                if( currentArcStateProcessors.size( ) == 0 )
                {
                    throw std::runtime_error( "Error when resetting multi-arc states, no state processors found" );
                }
                for( auto itr : currentArcStateProcessors )
                {
                    singleArcIntegratedStatesProcessors[ itr.first ].push_back( itr.second );
                }
            }

            std::map< IntegratedStateType, std::shared_ptr< MultiArcIntegratedStateProcessor< TimeType, StateScalarType > > >
                    multiArcStateProcessors = createMultiArcIntegratedStateProcessors(
                            bodies_, propagationResults_->getArcStartTimes( ), singleArcIntegratedStatesProcessors );
            for( auto itr : multiArcStateProcessors )
            {
                itr.second->setInterpolatorSettings( multiArcPropagatorSettings_->getOutputSettings( )->getInterpolatorSettings( ) );
                itr.second->processIntegratedMultiArcStates(
                        propagationResults_->getConcatenatedEquationsOfMotionResults(
                                multiArcPropagatorSettings_->getOutputSettings( )->getClearNumericalSolutions( ) ),
                        propagationResults_->getArcStartTimes( ) );
            }
        }
        catch( const std::exception& caughtException )
        {
            std::cerr << "Error occured when post-processing mulyi-arc integration results, and seting integrated states in "
                         "environment, caught error is: "
                      << std::endl
                      << std::endl;
            std::cerr << caughtException.what( ) << std::endl << std::endl;
            std::cerr << "The problem may be that there is an insufficient number of data points (epochs) at which propagation results "
                         "are produced for one or more arcs"
                      << std::endl;
            if( multiArcPropagatorSettings_->getOutputSettings( )->getClearNumericalSolutions( ) )
            {
                propagationResults_->clearSolutionMaps( );
            }
        }
    }
    else if( multiArcPropagatorSettings_->getOutputSettings( )->getClearNumericalSolutions( ) )
    {
        propagationResults_->clearSolutionMaps( );
    }

    // Reset dependent variables interface
    if( multiArcPropagatorSettings_->getOutputSettings( )->getUpdateDependentVariableInterpolator( ) )
    {
        propagationResults_->updateDependentVariableInterface( );
    }
}

template< typename StateScalarType, typename TimeType >
void MultiArcDynamicsSimulator< StateScalarType, TimeType >::printPrePropagationMessages( )
{
    if( multiArcPropagatorSettings_->getOutputSettings( )->printAnyOutput( ) )
    {
        std::cout << multiArcPropagatorSettings_->getOutputSettings( )->getPropagationStartHeader( ) << std::endl << std::endl;
    }
}

template< typename StateScalarType, typename TimeType >
void MultiArcDynamicsSimulator< StateScalarType, TimeType >::printPostPropagationMessages( )
{
    if( multiArcPropagatorSettings_->getOutputSettings( )->printAnyOutput( ) )
    {
        std::cout << multiArcPropagatorSettings_->getOutputSettings( )->getPropagationEndHeader( ) << std::endl << std::endl;
    }
}

template< typename StateScalarType, typename TimeType >
std::vector< double > MultiArcDynamicsSimulator< StateScalarType, TimeType >::getArcStartTimes( )
{
    return propagationResults_->getArcStartTimes( );
}

template< typename StateScalarType, typename TimeType >
std::vector< double > MultiArcDynamicsSimulator< StateScalarType, TimeType >::getArcEndTimes( )
{
    return propagationResults_->getArcEndTimes( );
}

template< typename StateScalarType, typename TimeType >
std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >
MultiArcDynamicsSimulator< StateScalarType, TimeType >::getEquationsOfMotionNumericalSolution( )
{
    return propagationResults_->getConcatenatedEquationsOfMotionResults( );
}

template< typename StateScalarType, typename TimeType >
std::vector< std::map< TimeType, Eigen::VectorXd > > MultiArcDynamicsSimulator< StateScalarType, TimeType >::getDependentVariableHistory( )
{
    return propagationResults_->getConcatenatedDependentVariableResults( );
}

template< typename StateScalarType, typename TimeType >
std::vector< std::map< TimeType, double > > MultiArcDynamicsSimulator< StateScalarType, TimeType >::getCumulativeComputationTimeHistory( )
{
    return propagationResults_->getConcatenatedCumulativeComputationTimeHistory( );
}

template< typename StateScalarType, typename TimeType >
std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >
MultiArcDynamicsSimulator< StateScalarType, TimeType >::getEquationsOfMotionNumericalSolutionBase( )
{
    return getEquationsOfMotionNumericalSolution( );
}

template< typename StateScalarType, typename TimeType >
std::vector< std::map< TimeType, Eigen::VectorXd > >
MultiArcDynamicsSimulator< StateScalarType, TimeType >::getDependentVariableNumericalSolutionBase( )
{
    return getDependentVariableHistory( );
}

template< typename StateScalarType, typename TimeType >
std::vector< std::map< TimeType, double > >
MultiArcDynamicsSimulator< StateScalarType, TimeType >::getCumulativeComputationTimeHistoryBase( )
{
    return getCumulativeComputationTimeHistory( );
}

template< typename StateScalarType, typename TimeType >
bool MultiArcDynamicsSimulator< StateScalarType, TimeType >::integrationCompletedSuccessfully( ) const
{
    return propagationResults_->integrationCompletedSuccessfully( );
}

template< typename StateScalarType, typename TimeType >
std::vector< std::shared_ptr< PropagationTerminationDetails > >
MultiArcDynamicsSimulator< StateScalarType, TimeType >::getPropagationTerminationReasons( )
{
    return propagationResults_->getConcatenatedTerminationReasons( );
}

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_MULTIARCDYNAMICSSIMULATORIMPLEMENTATION_H
