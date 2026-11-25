/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SINGLEARCDYNAMICSSIMULATORIMPLEMENTATION_H
#define TUDAT_SINGLEARCDYNAMICSSIMULATORIMPLEMENTATION_H

#include "tudat/simulation/propagation_setup/singleArcDynamicsSimulator.h"

namespace tudat
{

namespace propagators
{

template< typename StateScalarType, typename TimeType >
SingleArcDynamicsSimulator< StateScalarType, TimeType >::SingleArcDynamicsSimulator(
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
        const bool areEquationsOfMotionToBeIntegrated,
        const PredefinedSingleArcStateDerivativeModels< StateScalarType, TimeType >& predefinedStateDerivativeModels,
        const bool isPartOfMultiArc ):
    DynamicsSimulator< StateScalarType, TimeType >( bodies, propagatorSettings ),
    propagatorSettings_( std::dynamic_pointer_cast< SingleArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ) )
{
    // Check consistency of input settings
    if( propagatorSettings == nullptr )
    {
        throw std::runtime_error( "Error in dynamics simulator, propagator settings not defined." );
    }
    else if( propagatorSettings_ == nullptr )
    {
        throw std::runtime_error( "Error in dynamics simulator, input must be single-arc." );
    }
    else
    {
        // Retrieve output and integrator settings TODO: no need to set as member variables; can just retrieve from propagatorSettings_
        outputSettings_ = propagatorSettings_->getOutputSettingsWithCheck( );
        integratorSettings_ = propagatorSettings_->getIntegratorSettings( );
    }
    if( integratorSettings_ == nullptr )
    {
        throw std::runtime_error( "Error in dynamics simulator, integrator settings not defined." );
    }
    checkPropagatedStatesFeasibility( propagatorSettings_, bodies_, isPartOfMultiArc );

    // Create objects that reset the environment (e.g. ephemerides) after propagation is required
    if( propagatorSettings_->getOutputSettings( )->getCreateStateProcessors( ) )
    {
        createAndSetIntegratedStateProcessors( );
    }

    // Create object that updates the environment during propagation
    try
    {
        environmentUpdater_ = createEnvironmentUpdaterForDynamicalEquations< StateScalarType, TimeType >( propagatorSettings_, bodies_ );
    }
    catch( const std::runtime_error& error )
    {
        throw std::runtime_error( "Error when creating environment updater: " + std::string( error.what( ) ) );
    }

    // Create object that calculates the complete state derivatives
    if( predefinedStateDerivativeModels.stateDerivativeModels_.size( ) == 0 )
    {
        dynamicsStateDerivative_ = std::make_shared< DynamicsStateDerivativeModel< TimeType, StateScalarType > >(
                createStateDerivativeModels< StateScalarType, TimeType >(
                        propagatorSettings_, bodies_, propagatorSettings_->getInitialTime( ) ),
                std::bind( &EnvironmentUpdater< StateScalarType, TimeType >::updateEnvironment,
                           environmentUpdater_,
                           std::placeholders::_1,
                           std::placeholders::_2,
                           std::placeholders::_3 ) );
    }
    else
    {
        dynamicsStateDerivative_ = std::make_shared< DynamicsStateDerivativeModel< TimeType, StateScalarType > >(
                predefinedStateDerivativeModels.stateDerivativeModels_,
                std::bind( &EnvironmentUpdater< StateScalarType, TimeType >::updateEnvironment,
                           environmentUpdater_,
                           std::placeholders::_1,
                           std::placeholders::_2,
                           std::placeholders::_3 ) );
    }
    stateDerivativeFunction_ = std::bind( &DynamicsStateDerivativeModel< TimeType, StateScalarType >::computeStateDerivative,
                                          dynamicsStateDerivative_,
                                          std::placeholders::_1,
                                          std::placeholders::_2 );

    // Create object that determines if the propagation is to be terminated
    propagationTerminationCondition_ = createPropagationTerminationConditions( propagatorSettings_->getTerminationSettings( ),
                                                                               bodies_,
                                                                               integratorSettings_->initialTimeStep_,
                                                                               dynamicsStateDerivative_->getStateDerivativeModels( ),
                                                                               predefinedStateDerivativeModels.stateDerivativePartials_ );

    sequentialPropagation_ = true;
    if( propagationTerminationCondition_->getTerminationType( ) == non_sequential_stopping_condition )
    {
        sequentialPropagation_ = false;
        if( integratorSettings_->initialTimeStep_ < 0.0 )
        {
            throw std::runtime_error(
                    "Error when using non-sequential propagation, the initial integrator time step must be positive (first provided "
                    "for forward leg, "
                    "conversion to negative time step for backward leg is automatic)." );
        }
    }

    std::map< IntegratedStateType, std::vector< std::tuple< std::string, std::string, PropagatorType > > > integratedStateAndBodyList =
            getIntegratedTypeAndBodyList( propagatorSettings_ );

    std::map< std::pair< int, int >, std::string > dependentVariableIds_;
    std::map< std::pair< int, int >, std::shared_ptr< SingleDependentVariableSaveSettings > > orderedDependentVariableSettings_;

    // Create functions that compute the dependent variables
    if( propagatorSettings_->getDependentVariablesToSave( ).size( ) > 0 )
    {
        std::pair< std::function< Eigen::VectorXd( ) >, std::map< std::pair< int, int >, std::string > > dependentVariableData =
                createDependentVariableListFunction< TimeType, StateScalarType >( propagatorSettings_->getDependentVariablesToSave( ),
                                                                                  bodies_,
                                                                                  orderedDependentVariableSettings_,
                                                                                  dynamicsStateDerivative_->getStateDerivativeModels( ),
                                                                                  predefinedStateDerivativeModels.stateDerivativePartials_,
                                                                                  predefinedStateDerivativeModels.parametersToEstimate_ );
        dependentVariablesFunctions_ = dependentVariableData.first;
        dependentVariableIds_ = dependentVariableData.second;
    }

    // Create object that will contain and process the propagation results
    std::shared_ptr< SingleArcDependentVariablesInterface< TimeType > > dependentVariableInterface =
            std::make_shared< SingleArcDependentVariablesInterface< TimeType > >(
                    std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, Eigen::VectorXd > >( ),
                    propagatorSettings_->getDependentVariablesToSave( ),
                    dependentVariableIds_,
                    orderedDependentVariableSettings_,
                    bodies );

    propagationResults_ = std::make_shared< SingleArcSimulationResults< StateScalarType, TimeType > >(
            integratedStateAndBodyList,
            propagatorSettings_->getOutputSettingsWithCheck( ),
            std::bind( &DynamicsStateDerivativeModel< TimeType, StateScalarType >::convertNumericalStateSolutionsToOutputSolutions,
                       dynamicsStateDerivative_,
                       std::placeholders::_1,
                       std::placeholders::_2 ),
            dependentVariableInterface,
            sequentialPropagation_ );

    // Integrate equations of motion if required.
    if( areEquationsOfMotionToBeIntegrated )
    {
        integrateEquationsOfMotion( propagatorSettings_->getInitialStates( ) );
    }
}

template< typename StateScalarType, typename TimeType >
void SingleArcDynamicsSimulator< StateScalarType, TimeType >::integrateEquationsOfMotion(
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& initialStates )
{
    integrateEquationsOfMotion< SingleArcSimulationResults< StateScalarType, TimeType > >(
            dynamicsStateDerivative_->convertFromOutputSolution( initialStates, propagatorSettings_->getInitialTime( ) ),
            propagationResults_ );
}

template< typename StateScalarType, typename TimeType >
void SingleArcDynamicsSimulator< StateScalarType, TimeType >::integrate(
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& initialStates )
{
    integrateEquationsOfMotion( initialStates );
}

template< typename StateScalarType, typename TimeType >
template< typename SimulationResults >
void SingleArcDynamicsSimulator< StateScalarType, TimeType >::integrateEquationsOfMotion(
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& processedInitialState,
        const std::shared_ptr< SimulationResults > propagationResults )
{
    performPropagationPreProcessingSteps( propagationResults );
    propagateDynamics< SimulationResults >( processedInitialState, propagationResults );
    performPropagationPostProcessingSteps( propagationResults );
}

template< typename StateScalarType, typename TimeType >
std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >
SingleArcDynamicsSimulator< StateScalarType, TimeType >::getEquationsOfMotionNumericalSolutionBase( )
{
    return std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >(
            { getEquationsOfMotionNumericalSolution( ) } );
}

template< typename StateScalarType, typename TimeType >
std::vector< std::map< TimeType, Eigen::VectorXd > >
SingleArcDynamicsSimulator< StateScalarType, TimeType >::getDependentVariableNumericalSolutionBase( )
{
    return std::vector< std::map< TimeType, Eigen::VectorXd > >( { getDependentVariableHistory( ) } );
}

template< typename StateScalarType, typename TimeType >
std::vector< std::map< TimeType, double > >
SingleArcDynamicsSimulator< StateScalarType, TimeType >::getCumulativeComputationTimeHistoryBase( )
{
    return std::vector< std::map< TimeType, double > >( { getCumulativeComputationTimeHistory( ) } );
}

template< typename StateScalarType, typename TimeType >
double SingleArcDynamicsSimulator< StateScalarType, TimeType >::getInitialPropagationTime( )
{
    return propagatorSettings_->getInitialTime( );
}

template< typename StateScalarType, typename TimeType >
void SingleArcDynamicsSimulator< StateScalarType, TimeType >::resetInitialPropagationTime( const double initialPropagationTime )
{
    propagatorSettings_->resetInitialTime( initialPropagationTime );
}

template< typename StateScalarType, typename TimeType >
void SingleArcDynamicsSimulator< StateScalarType, TimeType >::resetPropagationTerminationConditions( )
{
    propagationTerminationCondition_ = createPropagationTerminationConditions( propagatorSettings_->getTerminationSettings( ),
                                                                               bodies_,
                                                                               integratorSettings_->initialTimeStep_,
                                                                               dynamicsStateDerivative_->getStateDerivativeModels( ) );
}

template< typename StateScalarType, typename TimeType >
void SingleArcDynamicsSimulator< StateScalarType, TimeType >::processNumericalEquationsOfMotionSolution( )
{
    if( outputSettings_->getSetIntegratedResult( ) )
    {
        try
        {
            for( auto processorIterator : integratedStateProcessors_ )
            {
                processorIterator.second->setInterpolatorSettings( outputSettings_->getInterpolatorSettings( ) );
            }
            // Create and set interpolators for ephemerides
            resetIntegratedStates( propagationResults_->equationsOfMotionNumericalSolution_, integratedStateProcessors_ );
        }
        catch( const std::exception& caughtException )
        {
            std::cerr << "Error occured when post-processing single-arc integration results, and seting integrated states in "
                         "environment, caught error is: "
                      << std::endl
                      << std::endl;
            std::cerr << caughtException.what( ) << std::endl << std::endl;
            std::cerr << "The problem may be that there is an insufficient number of data points (epochs) at which propagation results "
                         "are produced. Integrated results are given at" +
                            std::to_string( propagationResults_->equationsOfMotionNumericalSolution_.size( ) ) + " epochs"
                      << std::endl;
        }

        // Clear numerical solution if so required.
        if( propagatorSettings_->getOutputSettings( )->getClearNumericalSolutions( ) )
        {
            propagationResults_->clearSolutionMaps( );
        }

        for( auto bodyIterator : bodies_.getMap( ) )
        {
            bodyIterator.second->updateConstantEphemerisDependentMemberQuantities( );
        }
    }
    else if( propagatorSettings_->getOutputSettings( )->getClearNumericalSolutions( ) )
    {
        propagationResults_->clearSolutionMaps( );
    }

    if( propagatorSettings_->getOutputSettings( )->getUpdateDependentVariableInterpolator( ) )
    {
        propagationResults_->updateDependentVariableInterface( );
    }
}

template< typename StateScalarType, typename TimeType >
void SingleArcDynamicsSimulator< StateScalarType, TimeType >::suppressDependentVariableDataPrinting( )
{
    outputSettings_->getPrintSettings( )->setPrintDependentVariableData( false );
}

template< typename StateScalarType, typename TimeType >
void SingleArcDynamicsSimulator< StateScalarType, TimeType >::enableDependentVariableDataPrinting( )
{
    outputSettings_->getPrintSettings( )->setPrintDependentVariableData( true );
}

template< typename StateScalarType, typename TimeType >
void SingleArcDynamicsSimulator< StateScalarType, TimeType >::createAndSetIntegratedStateProcessors( )
{
    frameManager_ = simulation_setup::createFrameManager( bodies_.getMap( ) );
    integratedStateProcessors_ =
            createIntegratedStateProcessors< TimeType, StateScalarType >( propagatorSettings_, bodies_, frameManager_ );
}

template< typename StateScalarType, typename TimeType >
std::map< double, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >
SingleArcDynamicsSimulator< StateScalarType, TimeType >::getEquationsOfMotionNumericalSolutionDouble( )
{
    return propagationResults_->template getEquationsOfMotionNumericalSolutionTemplated< double >( );
}

template< typename StateScalarType, typename TimeType >
std::map< double, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >
SingleArcDynamicsSimulator< StateScalarType, TimeType >::getEquationsOfMotionNumericalSolutionRawDouble( )
{
    return propagationResults_->template getEquationsOfMotionNumericalSolutionRawTemplated< double >( );
}

template< typename StateScalarType, typename TimeType >
std::map< double, Eigen::VectorXd > SingleArcDynamicsSimulator< StateScalarType, TimeType >::getDependentVariableHistoryDouble( )
{
    return propagationResults_->template getDependentVariableHistoryTemplated< double >( );
}

template< typename StateScalarType, typename TimeType >
std::map< double, double > SingleArcDynamicsSimulator< StateScalarType, TimeType >::getCumulativeComputationTimeHistoryDouble( )
{
    return propagationResults_->getCumulativeComputationTimeHistory( );
}

template< typename StateScalarType, typename TimeType >
std::map< double, unsigned int > SingleArcDynamicsSimulator< StateScalarType, TimeType >::getCumulativeNumberOfFunctionEvaluationsDouble( )
{
    return propagationResults_->getCumulativeNumberOfFunctionEvaluations( );
}

template< typename StateScalarType, typename TimeType >
bool SingleArcDynamicsSimulator< StateScalarType, TimeType >::integrationCompletedSuccessfully( ) const
{
    return ( propagationResults_->propagationTerminationReason_->getPropagationTerminationReason( ) == termination_condition_reached );
}

template< typename StateScalarType, typename TimeType >
std::map< std::pair< int, int >, std::string > SingleArcDynamicsSimulator< StateScalarType, TimeType >::getDependentVariableIds( )
{
    return propagationResults_->getDependentVariableId( );
}

template< typename StateScalarType, typename TimeType >
template< typename InputTimeType, typename InputStateScalarType >
std::map< InputTimeType, Eigen::VectorXd >
SingleArcDynamicsSimulator< StateScalarType, TimeType >::evaluateDependentVariablesAlongTrajectory(
        const std::map< InputTimeType, Eigen::Matrix< InputStateScalarType, Eigen::Dynamic, 1 > >& stateHistory )
{
    std::map< InputTimeType, Eigen::VectorXd > dependentVariables_;

    // Integrate equations of motion numerically.
    simulation_setup::setAreBodiesInPropagation( bodies_, true );

    for( auto it : stateHistory )
    {
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > unprocessedState =
                dynamicsStateDerivative_->convertFromOutputSolution( it.second.template cast< StateScalarType >( ), it.first );
        dynamicsStateDerivative_->computeStateDerivative( static_cast< double >( it.first ), unprocessedState );
        dependentVariables_[ it.first ] = dependentVariablesFunctions_( );
    }

    simulation_setup::setAreBodiesInPropagation( bodies_, false );
    return dependentVariables_;
}

template< typename StateScalarType, typename TimeType >
template< typename SimulationResults >
void SingleArcDynamicsSimulator< StateScalarType, TimeType >::propagateDynamics(
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, SimulationResults::number_of_columns >& processedInitialState,
        const std::shared_ptr< SimulationResults > propagationResults )
{
    // Integrate equations of motion numerically.
    simulation_setup::setAreBodiesInPropagation( bodies_, true );
    dynamicsStateDerivative_->updateStateDerivativeModelSettings(
            processedInitialState.block( 0, processedInitialState.cols( ) - 1, processedInitialState.rows( ), 1 ) );

    if( sequentialPropagation_ )
    {
        propagateUntilTermination( processedInitialState, propagationResults, propagationTerminationCondition_ );
    }
    else
    {
        std::shared_ptr< NonSequentialPropagationTerminationCondition > nonSequentialTerminations =
                std::dynamic_pointer_cast< NonSequentialPropagationTerminationCondition >( propagationTerminationCondition_ );
        propagateUntilTermination(
                processedInitialState, propagationResults, nonSequentialTerminations->getForwardPropagationTerminationCondition( ) );

        integratorSettings_->initialTimeStep_ *= -1.0;
        propagateUntilTermination(
                processedInitialState, propagationResults, nonSequentialTerminations->getBackwardPropagationTerminationCondition( ) );
        integratorSettings_->initialTimeStep_ *= -1.0;
    }

    simulation_setup::setAreBodiesInPropagation( bodies_, false );
}

template< typename StateScalarType, typename TimeType >
template< typename SimulationResults >
void SingleArcDynamicsSimulator< StateScalarType, TimeType >::propagateUntilTermination(
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, SimulationResults::number_of_columns >& processedInitialState,
        const std::shared_ptr< SimulationResults > propagationResults,
        const std::shared_ptr< PropagationTerminationCondition > terminationCondition )
{
    SingleArcPropagationController< Eigen::Matrix< StateScalarType, Eigen::Dynamic, SimulationResults::number_of_columns >, TimeType >
            propagationController( dynamicsStateDerivative_,
                                   integratorSettings_,
                                   terminationCondition,
                                   dependentVariablesFunctions_,
                                   propagatorSettings_->getOutputSettings( ) );
    propagationController.initialize( processedInitialState, propagatorSettings_->getInitialTime( ) );
    propagationController.propagate( );
    propagationController.writeResults( propagationResults );
}

template< typename StateScalarType, typename TimeType >
template< typename SimulationResults >
void SingleArcDynamicsSimulator< StateScalarType, TimeType >::performPropagationPreProcessingSteps(
        const std::shared_ptr< SimulationResults > propagationResults )
{
    // Reset functions
    dynamicsStateDerivative_->setPropagationSettings( std::vector< IntegratedStateType >( ), true, SimulationResults::is_variational );
    dynamicsStateDerivative_->resetFunctionEvaluationCounter( );
    dynamicsStateDerivative_->resetCumulativeFunctionEvaluationCounter( );
    resetPropagationTerminationConditions( );

    // Empty solution maps
    propagationResults->reset( );

    PropagationPrintingInterface< SimulationResults, StateScalarType, TimeType >::printSingleArcPrePropagationMessages(
            outputSettings_->getPrintSettings( ), outputSettings_->getPropagationStartHeader( ), propagationResults );
}

template< typename StateScalarType, typename TimeType >
template< typename SimulationResults >
void SingleArcDynamicsSimulator< StateScalarType, TimeType >::performPropagationPostProcessingSteps(
        const std::shared_ptr< SimulationResults > propagationResults )
{
    // Retrieve number of cumulative function evaluations
    propagationResults->finalizePropagation( dynamicsStateDerivative_->getCumulativeNumberOfFunctionEvaluations( ) );
    PropagationPrintingInterface< SimulationResults, StateScalarType, TimeType >::printSingleArcPostPropagationMessages(
            outputSettings_->getPrintSettings( ), outputSettings_->getPropagationEndHeader( ), propagationResults );
    processNumericalEquationsOfMotionSolution( );
}

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_SINGLEARCDYNAMICSSIMULATORIMPLEMENTATION_H
