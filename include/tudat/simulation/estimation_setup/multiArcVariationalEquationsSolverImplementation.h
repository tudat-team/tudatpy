/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_MULTIARCVARIATIONALEQUATIONSSOLVERIMPLEMENTATION_H
#define TUDAT_MULTIARCVARIATIONALEQUATIONSSOLVERIMPLEMENTATION_H

#include "tudat/simulation/estimation_setup/multiArcVariationalEquationsSolver.h"

namespace tudat
{

namespace propagators
{

template< typename StateScalarType, typename TimeType >
MultiArcVariationalEquationsSolver< StateScalarType, TimeType >::MultiArcVariationalEquationsSolver(
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< MultiArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
        const bool integrateEquationsOnCreation ):
    VariationalEquationsSolver< StateScalarType, TimeType >(
            bodies,
            parametersToEstimate,
            propagatorSettings != nullptr ? propagatorSettings->getOutputSettingsWithCheck( )->getClearNumericalSolutions( ) : false ),
    propagatorSettings_( propagatorSettings ),
    resetMultiArcDynamicsAfterPropagation_(
            propagatorSettings != nullptr ? propagatorSettings->getOutputSettingsWithCheck( )->getSetIntegratedResult( ) : false )
{
    if( std::dynamic_pointer_cast< MultiArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ) == nullptr )
    {
        throw std::runtime_error( "Error when making multi-arc variational equations solver, input is single-arc" );
    }
    checkMultiArcPropagatorSettingsAndParameterEstimationConsistency( propagatorSettings_,
                                                                      parametersToEstimate,
                                                                      propagatorSettings->getArcStartTimes( ),
                                                                      estimatedBodiesPerArc_,
                                                                      arcIndicesPerBody_,
                                                                      areEstimatedBodiesDifferentPerArc_ );

    if( areEstimatedBodiesDifferentPerArc_ && resetMultiArcDynamicsAfterPropagation_ )
    {
        //            throw std::runtime_error( "Error in multi-arc variational equations solver, boolean
        //            resetMultiArcDynamicsAfterPropagation should be "
        //                                      "set to false when the bodies to be estimated differ from one arc to another." );
    }

    arcWiseParametersToEstimate_.clear( );
    estimatable_parameters::getParametersToEstimatePerArcTest( parametersToEstimate,
                                                               arcWiseParametersToEstimate_,
                                                               propagatorSettings->getArcStartTimes( ),
                                                               estimatedBodiesPerArc_,
                                                               arcIndicesPerBody_ );

    parameterVectorSize_ = 0;
    stateTransitionMatrixSize_ = 0;

    for( unsigned int arc = 0; arc < estimatedBodiesPerArc_.size( ); arc++ )
    {
        arcWiseStateTransitionMatrixSize_.push_back(
                estimatable_parameters::getSingleArcInitialDynamicalStateParameterSetSize( parametersToEstimate, arc ) );
        arcWiseParameterVectorSize_.push_back( estimatable_parameters::getSingleArcParameterSetSize( parametersToEstimate, arc ) );
    }

    dynamicsSimulator_ = std::make_shared< MultiArcDynamicsSimulator< StateScalarType, TimeType > >( bodies, propagatorSettings, false );

    std::vector< std::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > > singleArcDynamicsSimulators =
            dynamicsSimulator_->getSingleArcDynamicsSimulators( );

    std::vector< std::shared_ptr< SingleArcVariationalSimulationResults< StateScalarType, TimeType > > > singleArcVariationalResults;
    for( unsigned int i = 0; i < singleArcDynamicsSimulators.size( ); i++ )
    {
        dynamicsStateDerivatives_.push_back( singleArcDynamicsSimulators.at( i )->getDynamicsStateDerivative( ) );

        // Create variational equations objects.
        std::map< IntegratedStateType, orbit_determination::StateDerivativePartialsMap > stateDerivativePartials =
                simulation_setup::createStateDerivativePartials< StateScalarType, TimeType >(
                        dynamicsStateDerivatives_.at( i )->getStateDerivativeModels( ), bodies, arcWiseParametersToEstimate_[ i ] );

        std::shared_ptr< VariationalEquations > variationalEquationsObject_ = std::make_shared< VariationalEquations >(
                stateDerivativePartials,
                arcWiseParametersToEstimate_[ i ],
                dynamicsStateDerivatives_.at( i )->getStateTypeStartIndices( ),
                i,
                arcIndicesPerBody_[ i ],
                simulation_setup::getEphemerisOrigins( bodies ),
                getPropagationTimeInterval( propagatorSettings_->getSingleArcSettings( ).at( i ) ) );

        dynamicsStateDerivatives_.at( i )->addVariationalEquations( variationalEquationsObject_ );

        // Create object that will contain and process the propagation results, with the propagated number of columns
        singleArcVariationalResults.push_back( std::make_shared< SingleArcVariationalSimulationResults< StateScalarType, TimeType > >(
                singleArcDynamicsSimulators.at( i )->getSingleArcPropagationResults( ),
                arcWiseStateTransitionMatrixSize_.at( i ),
                arcWiseParameterVectorSize_.at( i ) - arcWiseStateTransitionMatrixSize_.at( i ) ) );
        singleArcVariationalResults.back( )->setReducedSensitivityMatrixColumns(
                variationalEquationsObject_->getReducedSensitivityMatrixColumns( ) );
        arcWiseParameterVectorSize_.at( i ) = variationalEquationsObject_->getNumberOfParameterValues( );
    }
    variationalPropagationResults_ =
            std::make_shared< MultiArcSimulationResults< SingleArcVariationalSimulationResults, StateScalarType, TimeType > >(
                    singleArcVariationalResults,
                    dynamicsSimulator_->getMultiArcPropagationResults( )->getMultiArcDependentVariablesInterface( ) );

    numberOfArcs_ = dynamicsStateDerivatives_.size( );
    // Integrate variational equations from initial state estimate.
    if( integrateEquationsOnCreation )
    {
        integrateVariationalAndDynamicalEquations( propagatorSettings_->getInitialStateList( ), 1 );
    }
    else
    {
        stateTransitionInterface_ = std::make_shared< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface< StateScalarType > >(
                propagatorSettings_->getArcStartTimes( ),
                parametersToEstimate_,
                propagatorSettings_->getSingleArcSettings( ).at( 0 )->getConventionalStateSize( ),
                parametersToEstimate_->getParameterSetSize( ) );
        setReducedSensitivityMatrixColumnsInInterface( );
    }
}

template< typename StateScalarType, typename TimeType >
void MultiArcVariationalEquationsSolver< StateScalarType, TimeType >::integrateDynamicalEquationsOfMotionOnly(
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialStateEstimate )
{
    dynamicsSimulator_->integrateEquationsOfMotion( initialStateEstimate );
}

template< typename StateScalarType, typename TimeType >
void MultiArcVariationalEquationsSolver< StateScalarType, TimeType >::integrateDynamicalEquationsOfMotionOnly(
        const std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& initialStateEstimate )
{
    dynamicsSimulator_->integrateEquationsOfMotion( initialStateEstimate );
}

template< typename StateScalarType, typename TimeType >
void MultiArcVariationalEquationsSolver< StateScalarType, TimeType >::integrateVariationalAndDynamicalEquations(
        const VectorType& concatenatedInitialStates,
        const bool integrateEquationsConcurrently )
{
    std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > splitInitialState;

    int currentIndex = 0;
    for( unsigned int i = 0; i < dynamicsSimulator_->getSingleArcDynamicsSimulators( ).size( ); i++ )
    {
        int currentSize =
                dynamicsSimulator_->getSingleArcDynamicsSimulators( ).at( i )->getPropagatorSettings( )->getConventionalStateSize( );
        splitInitialState.push_back( concatenatedInitialStates.block( currentIndex, 0, currentSize, 1 ) );
        currentIndex += currentSize;
    }

    if( currentIndex != concatenatedInitialStates.rows( ) )
    {
        throw std::runtime_error(
                "Error when doing multi-arc variational equation integration, "
                "input state vector size is incompatible with settings." );
    }
    integrateVariationalAndDynamicalEquations( splitInitialState, integrateEquationsConcurrently );
}

template< typename StateScalarType, typename TimeType >
std::shared_ptr< MultiArcInitialStateProvider< StateScalarType > >
MultiArcVariationalEquationsSolver< StateScalarType, TimeType >::getInitialStateProvider(
        const std::vector< VectorType >& initialStateEstimate )
{
    std::vector< std::pair< int, int > > variationalEquationsSize =
            utilities::mergeVectorsIntoVectorOfPairs( arcWiseStateTransitionMatrixSize_, arcWiseParameterVectorSize_ );
    return std::make_shared< MultiArcInitialStateProvider< StateScalarType > >( initialStateEstimate, variationalEquationsSize );
}

template< typename StateScalarType, typename TimeType >
void MultiArcVariationalEquationsSolver< StateScalarType, TimeType >::integrateVariationalAndDynamicalEquations(
        const std::vector< VectorType >& initialStateEstimate,
        const bool integrateEquationsConcurrently )
{
    // Propagate variational equations and equations of motion concurrently
    if( integrateEquationsConcurrently )
    {
        for( unsigned int i = 0; i < dynamicsStateDerivatives_.size( ); i++ )
        {
            dynamicsStateDerivatives_.at( i )->getVariationalEquationsCalculator( )->checkPropagationTimeInterval(
                    getPropagationTimeInterval( propagatorSettings_->getSingleArcSettings( ).at( i ) ) );
        }

        // Propagate dynamics and variational equations and store results in variationalPropagationResults_ object
        dynamicsSimulator_->template integrateEquationsOfMotion< MultiArcVariationalResults >(
                variationalPropagationResults_, getInitialStateProvider( initialStateEstimate ) );

        // Ensure consistency between parameters and propagator settings
        setPropagatorSettingsMultiArcStatesInEstimatedDynamicalParameters< StateScalarType, TimeType >( parametersToEstimate_,
                                                                                                        propagatorSettings_ );
    }

    // Reset solution for state transition and sensitivity matrices.
    if( propagatorSettings_->getOutputSettings( )->getSetIntegratedVariationalResult( ) )
    {
        resetVariationalEquationsInterpolators( );
    }
}

template< typename StateScalarType, typename TimeType >
std::shared_ptr< DynamicsSimulator< StateScalarType, TimeType > >
MultiArcVariationalEquationsSolver< StateScalarType, TimeType >::getDynamicsSimulatorBase( )
{
    return getDynamicsSimulator( );
}

template< typename StateScalarType, typename TimeType >
void MultiArcVariationalEquationsSolver< StateScalarType, TimeType >::resetParameterEstimate(
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > newParameterEstimate,
        const bool areVariationalEquationsToBeIntegrated )
{
    // Reset values of parameters.
    parametersToEstimate_->template resetParameterValues< StateScalarType >( newParameterEstimate );
    simulation_setup::setInitialStateVectorFromParameterSet< StateScalarType, TimeType >( parametersToEstimate_, propagatorSettings_ );

    for( int i = 0; i < numberOfArcs_; i++ )
    {
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > newParametersValues =
                propagatorSettings_->getSingleArcSettings( ).at( i )->getInitialStates( );
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > arcWiseParametersValues =
                arcWiseParametersToEstimate_.at( i )->template getFullParameterValues< StateScalarType >( );
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > newArcWiseParametersValues = arcWiseParametersValues;
        newArcWiseParametersValues.segment( 0, newParametersValues.size( ) ) = newParametersValues;
        arcWiseParametersToEstimate_.at( i )->resetParameterValues( newArcWiseParametersValues );
    }

    // Check if re-integration of variational equations is requested
    if( areVariationalEquationsToBeIntegrated )
    {
        // Integrate variational and state equations.
        this->integrateVariationalAndDynamicalEquations( propagatorSettings_->getInitialStates( ), 1 );
    }
    else
    {
        this->integrateDynamicalEquationsOfMotionOnly( propagatorSettings_->getInitialStates( ) );
    }
}

template< typename StateScalarType, typename TimeType >
std::vector< std::vector< std::map< double, Eigen::MatrixXd > > >
MultiArcVariationalEquationsSolver< StateScalarType, TimeType >::getNumericalVariationalEquationsSolution( )
{
    std::vector< std::vector< std::map< double, Eigen::MatrixXd > > > fullSolution;
    fullSolution.resize( numberOfArcs_ );
    for( int i = 0; i < numberOfArcs_; i++ )
    {
        fullSolution[ i ].push_back( variationalPropagationResults_->getSingleArcResults( ).at( i )->getStateTransitionSolution( ) );
        fullSolution[ i ].push_back( variationalPropagationResults_->getSingleArcResults( ).at( i )->getSensitivitySolution( ) );
    }
    std::cerr << "Warning, use of deprecated multi-arc getNumericalVariationalEquationsSolution is not recommended" << std::endl;
    return fullSolution;
}

template< typename StateScalarType, typename TimeType >
std::shared_ptr< SimulationResults< StateScalarType, TimeType > >
MultiArcVariationalEquationsSolver< StateScalarType, TimeType >::getVariationalPropagationResults( )
{
    return getMultiArcVariationalPropagationResults( );
}

template< typename StateScalarType, typename TimeType >
void MultiArcVariationalEquationsSolver< StateScalarType, TimeType >::resetVariationalEquationsInterpolators( )
{
    using namespace interpolators;

    // Allocate interpolator vectors
    std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
            stateTransitionMatrixInterpolators;
    std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > > sensitivityMatrixInterpolators;
    stateTransitionMatrixInterpolators.resize( variationalPropagationResults_->getSingleArcResults( ).size( ) );
    sensitivityMatrixInterpolators.resize( variationalPropagationResults_->getSingleArcResults( ).size( ) );

    // Create interpolators.
    std::vector< double > arcStartTimesToUse;
    std::vector< double > arcEndTimesToUse;

    for( unsigned int i = 0; i < variationalPropagationResults_->getSingleArcResults( ).size( ); i++ )
    {
        arcStartTimesToUse.push_back( dynamicsSimulator_->getArcStartTimes( ).at( i ) );
        arcEndTimesToUse.push_back( dynamicsSimulator_->getArcEndTimes( ).at( i ) );

        try
        {
            createStateTransitionAndSensitivityMatrixInterpolator(
                    stateTransitionMatrixInterpolators[ i ],
                    sensitivityMatrixInterpolators[ i ],
                    variationalPropagationResults_->getSingleArcResults( ).at( i )->getStateTransitionSolution( ),
                    variationalPropagationResults_->getSingleArcResults( ).at( i )->getReducedSensitivitySolution( ),
                    this->clearNumericalSolution_ );
        }
        catch( const std::exception& caughtException )
        {
            std::cerr << "Error occured when post-processing multi-arc variational equation integration results, and creating "
                         "interpolators in arc" +
                            std::to_string( i ) + ", caught error is: "
                      << std::endl
                      << std::endl;
            std::cerr << caughtException.what( ) << std::endl << std::endl;
            std::cerr << "The problem may be that there is an insufficient number of data points (epochs) at which propagation results "
                         "are produced for one or more arcs. Integrated results are given at" +
                            std::to_string( variationalPropagationResults_->getSingleArcResults( )
                                                    .at( 0 )
                                                    ->getStateTransitionSolution( )
                                                    .size( ) ) +
                            " epochs"
                      << std::endl;
        }
    }

    // Create stare transition matrix interface if needed, reset otherwise.
    if( stateTransitionInterface_ == nullptr )
    {
        stateTransitionInterface_ = std::make_shared< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface< StateScalarType > >(
                stateTransitionMatrixInterpolators,
                sensitivityMatrixInterpolators,
                propagatorSettings_->getArcStartTimes( ),
                arcStartTimesToUse,
                arcEndTimesToUse,
                parametersToEstimate_,
                propagatorSettings_->getSingleArcSettings( ).at( 0 )->getConventionalStateSize( ),
                parametersToEstimate_->getParameterSetSize( ),
                getArcWiseStatePartialAdditionIndices( ) );
        setReducedSensitivityMatrixColumnsInInterface( );
    }
    else
    {
        std::dynamic_pointer_cast< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface< StateScalarType > >(
                stateTransitionInterface_ )
                ->updateMatrixInterpolators( stateTransitionMatrixInterpolators,
                                             sensitivityMatrixInterpolators,
                                             arcStartTimesToUse,
                                             arcEndTimesToUse,
                                             getArcWiseStatePartialAdditionIndices( ) );
    }
}

template< typename StateScalarType, typename TimeType >
void MultiArcVariationalEquationsSolver< StateScalarType, TimeType >::setReducedSensitivityMatrixColumnsInInterface( )
{
    std::vector< std::vector< int > > reducedSensitivityMatrixColumns;
    for( unsigned int i = 0; i < dynamicsStateDerivatives_.size( ); i++ )
    {
        reducedSensitivityMatrixColumns.push_back(
                dynamicsStateDerivatives_.at( i )->getVariationalEquationsCalculator( )->getReducedSensitivityMatrixColumns( ) );
    }
    std::dynamic_pointer_cast< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface< StateScalarType > >(
            stateTransitionInterface_ )
            ->setReducedSensitivityMatrixColumns( reducedSensitivityMatrixColumns );
}

template< typename StateScalarType, typename TimeType >
std::vector< std::vector< std::pair< int, int > > >
MultiArcVariationalEquationsSolver< StateScalarType, TimeType >::getArcWiseStatePartialAdditionIndices( )
{
    std::vector< std::vector< std::pair< int, int > > > partialIndices;
    for( unsigned int i = 0; i < dynamicsStateDerivatives_.size( ); i++ )
    {
        partialIndices.push_back(
                dynamicsStateDerivatives_.at( i )->getVariationalEquationsCalculator( )->getStatePartialAdditionIndices( ) );
    }
    return partialIndices;
}

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_MULTIARCVARIATIONALEQUATIONSSOLVERIMPLEMENTATION_H
