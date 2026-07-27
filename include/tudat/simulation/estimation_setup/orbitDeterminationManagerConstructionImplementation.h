/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ORBITDETERMINATIONMANAGERCONSTRUCTIONIMPLEMENTATION_H
#define TUDAT_ORBITDETERMINATIONMANAGERCONSTRUCTIONIMPLEMENTATION_H

#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"
#include "tudat/basics/utilities.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/estimation_setup/createObservationManager.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"

namespace tudat
{

namespace simulation_setup
{

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
OrbitDeterminationManager< ObservationScalarType, TimeType, Dummy >::OrbitDeterminationManager(
        const SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > parametersToEstimate,
        const std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >& observationSettingsList,
        const std::shared_ptr< propagators::PropagatorSettings< ObservationScalarType > > propagatorSettings,
        const bool propagateOnCreation ):
    parametersToEstimate_( parametersToEstimate ), considerParameters_( parametersToEstimate_->getConsiderParameters( ) ), bodies_( bodies )
{
    initializeOrbitDeterminationManager( bodies, observationSettingsList, propagatorSettings, propagateOnCreation );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void OrbitDeterminationManager< ObservationScalarType, TimeType, Dummy >::initializeOrbitDeterminationManager(
        const SystemOfBodies& bodies,
        const std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >& observationSettingsList,
        const std::shared_ptr< propagators::PropagatorSettings< ObservationScalarType > > propagatorSettings,
        const bool propagateOnCreation )
{
    propagators::toggleIntegratedResultSettings< ObservationScalarType, TimeType >( propagatorSettings );

    // Detect whether consider parameters are included
    considerParametersIncluded_ = false;
    if( considerParameters_ != nullptr )
    {
        considerParametersIncluded_ = true;
    }

    // Retrieve size of estimated and consider parameters
    totalNumberParameters_ = parametersToEstimate_->getParameterSetSize( );
    numberEstimatedParameters_ = parametersToEstimate_->getEstimatedParameterSetSize( );
    numberConsiderParameters_ = 0;
    if( considerParameters_ != nullptr )
    {
        numberConsiderParameters_ = considerParameters_->getParameterSetSize( );
        if( !( ( totalNumberParameters_ - numberEstimatedParameters_ ) == numberConsiderParameters_ ) )
        {
            throw std::runtime_error( "Error when creating Estimator, number of consider parameters is inconsistent" );
        }
    }

    // Check if any dynamics is to be estimated
    std::map< propagators::IntegratedStateType, std::vector< std::pair< std::string, std::string > > > initialDynamicalStates =
            estimatable_parameters::getListOfInitialDynamicalStateParametersEstimate< ObservationScalarType >( parametersToEstimate_ );
    integrateAndEstimateOrbit_ = ( initialDynamicalStates.size( ) > 0 );

    propagatorSettings->getOutputSettingsBase( )->setUpdateDependentVariableInterpolator( true );
    if( integrateAndEstimateOrbit_ )
    {
        variationalEquationsSolver_ = simulation_setup::createVariationalEquationsSolver< ObservationScalarType, TimeType >(
                bodies, propagatorSettings, parametersToEstimate_, propagateOnCreation );
    }

    if( integrateAndEstimateOrbit_ )
    {
        stateTransitionAndSensitivityMatrixInterface_ = variationalEquationsSolver_->getStateTransitionMatrixInterface( );
    }
    else if( propagatorSettings == nullptr )
    {
        stateTransitionAndSensitivityMatrixInterface_ =
                createStateTransitionAndSensitivityMatrixInterface< ObservationScalarType, TimeType >(
                        propagatorSettings, parametersToEstimate_, 0, totalNumberParameters_ );
    }
    else
    {
        throw std::runtime_error( "Error, cannot parse propagator settings without estimating dynamics in OrbitDeterminationManager" );
    }

    // TODO correct this when moving dependent variable interface into results object
    if( std::dynamic_pointer_cast< propagators::HybridArcVariationalEquationsSolver< ObservationScalarType, TimeType > >(
                variationalEquationsSolver_ ) == nullptr )
    {
        dependentVariablesInterface_ = variationalEquationsSolver_->getDynamicsSimulatorBase( )->getDependentVariablesInterface( );
    }

    // Iterate over all observables and create observation managers.
    observationManagers_ = createObservationManagersBase( observationSettingsList,
                                                          bodies,
                                                          parametersToEstimate_,
                                                          stateTransitionAndSensitivityMatrixInterface_,
                                                          dependentVariablesInterface_ );

    // Set current parameter estimate from body initial states and parameter set.
    currentParameterEstimate_ = parametersToEstimate_->template getFullParameterValues< ObservationScalarType >( false );
    if( considerParametersIncluded_ )
    {
        considerParametersValues_ = considerParameters_->template getFullParameterValues< ObservationScalarType >( );
    }
    else
    {
        considerParametersValues_ = ParameterVectorType::Zero( 0 );
    }
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_ORBITDETERMINATIONMANAGERCONSTRUCTIONIMPLEMENTATION_H
