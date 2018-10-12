/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/SimulationSetup/PropagationSetup/createStateDerivativeModel.h"

namespace tudat
{

namespace propagators
{

//! Function to create an integrator to propagate the dynamics (in normalized units) in CR3BP
std::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::Vector6d > > createCR3BPIntegrator(
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const double massParameter,
        const Eigen::Vector6d& initialState )
{
    // Create state derivative model
    std::shared_ptr< StateDerivativeCircularRestrictedThreeBodyProblem > stateDerivativeModel =
            std::make_shared< StateDerivativeCircularRestrictedThreeBodyProblem  >( massParameter );
    std::function< Eigen::Vector6d( const double, const Eigen::Vector6d& ) > stateDerivativeFunction =
            std::bind( &StateDerivativeCircularRestrictedThreeBodyProblem::computeStateDerivative, stateDerivativeModel,
                       std::placeholders::_1, std::placeholders::_2 );

    // Create integrator object
    return numerical_integrators::createIntegrator< double, Eigen::Vector6d >(
                stateDerivativeFunction, initialState, integratorSettings );
}

//! Function to propagate the dynamics (in normalized units) in CR3BP
std::map< double, Eigen::Vector6d > performCR3BPIntegration(
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const double massParameter,
        const Eigen::Vector6d& initialState,
        const double finalTime )
{
    // Create integrator object
    std::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::Vector6d > > integrator =
            createCR3BPIntegrator( integratorSettings, massParameter, initialState );

    // Initialize return data map
    std::map< double, Eigen::Vector6d > stateHistory;

    // Store initial state and time
    double currentTime = integratorSettings->initialTime_;
    Eigen::Vector6d currentState = initialState;
    stateHistory[ currentTime ] = currentState;

    // Integrate to final time
    double timeStep = integratorSettings->initialTimeStep_;
    while( currentTime <= finalTime )
    {
        currentState = integrator->performIntegrationStep( timeStep );
        currentTime = integrator->getCurrentIndependentVariable( );
        timeStep = integrator->getNextStepSize( );
        stateHistory[ currentTime ] = currentState;
    }

    return stateHistory;

}

template std::vector< std::shared_ptr< SingleStateTypeDerivative< double, double > > > createStateDerivativeModels< double, double >(
        const std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const double propagationStartTime );
template std::shared_ptr< SingleStateTypeDerivative< double, double > > createStateDerivativeModel< double, double >(
        const std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const double propagationStartTime );

#if( BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS )
template std::vector< std::shared_ptr< SingleStateTypeDerivative< long double, double > > > createStateDerivativeModels< long double, double >(
        const std::shared_ptr< SingleArcPropagatorSettings< long double > > propagatorSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const double propagationStartTime );
template std::vector< std::shared_ptr< SingleStateTypeDerivative< double, Time > > > createStateDerivativeModels< double, Time >(
        const std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const Time propagationStartTime );
template std::vector< std::shared_ptr< SingleStateTypeDerivative< long double, Time > > > createStateDerivativeModels< long double, Time >(
        const std::shared_ptr< SingleArcPropagatorSettings< long double > > propagatorSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const Time propagationStartTime );
template std::shared_ptr< SingleStateTypeDerivative< long double, double > > createStateDerivativeModel< long double, double >(
        const std::shared_ptr< SingleArcPropagatorSettings< long double > > propagatorSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const double propagationStartTime );
template std::shared_ptr< SingleStateTypeDerivative< double, Time > > createStateDerivativeModel< double, Time >(
        const std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const Time propagationStartTime );
template std::shared_ptr< SingleStateTypeDerivative< long double, Time > > createStateDerivativeModel< long double, Time >(
        const std::shared_ptr< SingleArcPropagatorSettings< long double > > propagatorSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const Time propagationStartTime );
#endif

}

}
