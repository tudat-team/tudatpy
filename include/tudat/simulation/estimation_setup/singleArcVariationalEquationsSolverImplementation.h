/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SINGLEARCVARIATIONALEQUATIONSSOLVERIMPLEMENTATION_H
#define TUDAT_SINGLEARCVARIATIONALEQUATIONSSOLVERIMPLEMENTATION_H

#include "tudat/simulation/estimation_setup/singleArcVariationalEquationsSolver.h"

namespace tudat
{

namespace propagators
{

template< typename StateScalarType, typename TimeType >
SingleArcVariationalEquationsSolver< StateScalarType, TimeType >::SingleArcVariationalEquationsSolver(
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const bool integrateEquationsOnCreation ):
    VariationalEquationsSolver< StateScalarType, TimeType >(
            bodies,
            parametersToEstimate,
            propagatorSettings != nullptr && propagatorSettings->getOutputSettingsBase( ) != nullptr
                    ? propagatorSettings->getOutputSettingsBase( )->getClearNumericalSolutions( )
                    : false ),
    propagatorSettings_( std::dynamic_pointer_cast< SingleArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ) )
{
    // Check input consistency
    if( propagatorSettings_ == nullptr )
    {
        throw std::runtime_error( "Error in variational equations solver, input must be single-arc." );
    }
    else if( !checkPropagatorSettingsAndParameterEstimationConsistency< StateScalarType, TimeType >( propagatorSettings_,
                                                                                                     parametersToEstimate ) )
    {
        throw std::runtime_error(
                "Error when making single arc variational equations solver, estimated and propagated bodies are inconsistent." );
    }

    // Create state derivative models
    std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > stateDerivativeModels =
            createStateDerivativeModels( propagatorSettings_, bodies, propagatorSettings_->getInitialTime( ) );

    // Create state derivative partials
    std::map< IntegratedStateType, orbit_determination::StateDerivativePartialsMap > stateDerivativePartials =
            simulation_setup::createStateDerivativePartials< StateScalarType, TimeType >(
                    getStateDerivativeModelMapFromVector( stateDerivativeModels ), bodies, parametersToEstimate );

    // Create object that propagates the dynamics
    dynamicsSimulator_ = std::make_shared< SingleArcDynamicsSimulator< StateScalarType, TimeType > >(
            bodies,
            propagatorSettings_,
            false,
            PredefinedSingleArcStateDerivativeModels< StateScalarType, TimeType >(
                    stateDerivativeModels, stateDerivativePartials, parametersToEstimate_ ) );

    // Create variational equations evaluation objects.
    variationalEquationsObject_ =
            std::make_shared< VariationalEquations >( stateDerivativePartials,
                                                      parametersToEstimate_,
                                                      dynamicsSimulator_->getDynamicsStateDerivative( )->getStateTypeStartIndices( ),
                                                      -1,
                                                      std::map< std::string, int >( ),
                                                      simulation_setup::getEphemerisOrigins( bodies ),
                                                      getPropagationTimeInterval( propagatorSettings_ ) );
    dynamicsSimulator_->getDynamicsStateDerivative( )->addVariationalEquations( variationalEquationsObject_ );
    this->propagatedParameterVectorSize_ = variationalEquationsObject_->getNumberOfParameterValues( );

    // Create object that will contain and process the propagation results
    variationalPropagationResults_ = std::make_shared< SingleArcVariationalSimulationResults< StateScalarType, TimeType > >(
            dynamicsSimulator_->getSingleArcPropagationResults( ),
            this->stateTransitionMatrixSize_,
            this->parameterVectorSize_ - this->stateTransitionMatrixSize_ );
    variationalPropagationResults_->setReducedSensitivityMatrixColumns(
            variationalEquationsObject_->getReducedSensitivityMatrixColumns( ) );

    // Integrate variational equations from initial state estimate.
    if( integrateEquationsOnCreation )
    {
        if( integrateDynamicalAndVariationalEquationsConcurrently )
        {
            integrateVariationalAndDynamicalEquations( propagatorSettings_->getInitialStates( ), true );
        }
        else
        {
            integrateVariationalAndDynamicalEquations( propagatorSettings_->getInitialStates( ), false );
        }
    }
    else
    {
        createStateTransitionInterface( nullptr, nullptr );
    }
}

template< typename StateScalarType, typename TimeType >
void SingleArcVariationalEquationsSolver< StateScalarType, TimeType >::integrateDynamicalEquationsOfMotionOnly(
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialStateEstimate )
{
    dynamicsSimulator_->integrateEquationsOfMotion( initialStateEstimate );
}

template< typename StateScalarType, typename TimeType >
void SingleArcVariationalEquationsSolver< StateScalarType, TimeType >::integrateVariationalAndDynamicalEquations(
        const VectorType& initialStateEstimate,
        const bool integrateEquationsConcurrently )
{
    if( integrateEquationsConcurrently )
    {
        variationalEquationsObject_->checkPropagationTimeInterval( getPropagationTimeInterval( propagatorSettings_ ) );

        // Create initial conditions from new estimate.
        MatrixType initialVariationalState =
                this->createInitialConditions( dynamicsSimulator_->getDynamicsStateDerivative( )->convertFromOutputSolution(
                        initialStateEstimate, propagatorSettings_->getInitialTime( ) ) );

        // Propagate dynamics and variational equations
        dynamicsSimulator_->integrateEquationsOfMotion( initialVariationalState, variationalPropagationResults_ );
    }

    // Reset solution for state transition and sensitivity matrices.
    if( propagatorSettings_->getOutputSettings( )->getSetIntegratedVariationalResult( ) )
    {
        resetVariationalEquationsInterpolators( );
    }
}

template< typename StateScalarType, typename TimeType >
std::vector< std::map< double, Eigen::MatrixXd > >
SingleArcVariationalEquationsSolver< StateScalarType, TimeType >::getNumericalVariationalEquationsSolution( )
{
    std::cerr << "Warning, use of deprecated single-arc getNumericalVariationalEquationsSolution is not recommended" << std::endl;
    return std::vector< std::map< double, Eigen::MatrixXd > >( { getStateTransitionMatrixSolution( ), getSensitivityMatrixSolution( ) } );
}

template< typename StateScalarType, typename TimeType >
std::map< double, Eigen::MatrixXd >& SingleArcVariationalEquationsSolver< StateScalarType, TimeType >::getStateTransitionMatrixSolution( )
{
    return variationalPropagationResults_->getStateTransitionSolution( );
}

template< typename StateScalarType, typename TimeType >
std::map< double, Eigen::MatrixXd > SingleArcVariationalEquationsSolver< StateScalarType, TimeType >::getSensitivityMatrixSolution( )
{
    return variationalPropagationResults_->getSensitivitySolution( );
}

template< typename StateScalarType, typename TimeType >
const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >&
SingleArcVariationalEquationsSolver< StateScalarType, TimeType >::getEquationsOfMotionSolution( )
{
    return dynamicsSimulator_->getEquationsOfMotionNumericalSolution( );
}

template< typename StateScalarType, typename TimeType >
std::map< double, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >
SingleArcVariationalEquationsSolver< StateScalarType, TimeType >::getEquationsOfMotionSolutionDouble( )
{
    return dynamicsSimulator_->getEquationsOfMotionNumericalSolutionDouble( );
}

template< typename StateScalarType, typename TimeType >
std::shared_ptr< DynamicsSimulator< StateScalarType, TimeType > >
SingleArcVariationalEquationsSolver< StateScalarType, TimeType >::getDynamicsSimulatorBase( )
{
    return getDynamicsSimulator( );
}

template< typename StateScalarType, typename TimeType >
void SingleArcVariationalEquationsSolver< StateScalarType, TimeType >::resetParameterEstimate(
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > newParameterEstimate,
        const bool areVariationalEquationsToBeIntegrated )
{
    // Reset values of parameters.
    parametersToEstimate_->template resetParameterValues< StateScalarType >( newParameterEstimate );
    simulation_setup::setInitialStateVectorFromParameterSet< StateScalarType, TimeType >( parametersToEstimate_, propagatorSettings_ );

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
std::shared_ptr< SimulationResults< StateScalarType, TimeType > >
SingleArcVariationalEquationsSolver< StateScalarType, TimeType >::getVariationalPropagationResults( )
{
    return getSingleArcVariationalPropagationResults( );
}

template< typename StateScalarType, typename TimeType >
void SingleArcVariationalEquationsSolver< StateScalarType, TimeType >::resetVariationalEquationsInterpolators( )
{
    using namespace interpolators;
    using namespace utilities;

    // Create interpolators.
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > stateTransitionMatrixInterpolator;
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > sensitivityMatrixInterpolator;

    try
    {
        createStateTransitionAndSensitivityMatrixInterpolator( stateTransitionMatrixInterpolator,
                                                               sensitivityMatrixInterpolator,
                                                               variationalPropagationResults_->getStateTransitionSolution( ),
                                                               variationalPropagationResults_->getReducedSensitivitySolution( ),
                                                               this->clearNumericalSolution_ );
    }
    catch( const std::exception& caughtException )
    {
        std::cerr << "Error occured when post-processing single-arc variational equation integration results, and creating "
                     "interpolators, caught error is: "
                  << std::endl
                  << std::endl;
        std::cerr << caughtException.what( ) << std::endl << std::endl;
        std::cerr << "The problem may be that there is an insufficient number of data points (epochs) at which propagation results are "
                     "produced for one or more arcs. Integrated results are given at" +
                        std::to_string( variationalPropagationResults_->getStateTransitionSolution( ).size( ) ) + " epochs"
                  << std::endl;
    }

    // Create (if non-existent) or reset state transition matrix interface
    if( stateTransitionInterface_ == nullptr )
    {
        createStateTransitionInterface( stateTransitionMatrixInterpolator, sensitivityMatrixInterpolator );
    }
    else
    {
        std::dynamic_pointer_cast< SingleArcCombinedStateTransitionAndSensitivityMatrixInterface >( stateTransitionInterface_ )
                ->updateMatrixInterpolators( stateTransitionMatrixInterpolator,
                                             sensitivityMatrixInterpolator,
                                             variationalEquationsObject_->getStatePartialAdditionIndices( ) );
    }
}

template< typename StateScalarType, typename TimeType >
void SingleArcVariationalEquationsSolver< StateScalarType, TimeType >::createStateTransitionInterface(
        const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > stateTransitionMatrixInterpolator,
        const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > sensitivityMatrixInterpolator )
{
    std::shared_ptr< SingleArcCombinedStateTransitionAndSensitivityMatrixInterface > singleArcInterface =
            std::make_shared< SingleArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                    stateTransitionMatrixInterpolator,
                    sensitivityMatrixInterpolator,
                    propagatorSettings_->getConventionalStateSize( ),
                    parameterVectorSize_,
                    variationalEquationsObject_->getStatePartialAdditionIndices( ) );
    singleArcInterface->setReducedSensitivityMatrixColumns( variationalEquationsObject_->getReducedSensitivityMatrixColumns( ) );
    stateTransitionInterface_ = singleArcInterface;
}

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_SINGLEARCVARIATIONALEQUATIONSSOLVERIMPLEMENTATION_H
