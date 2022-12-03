/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/propagation_setup/createStateDerivativeModel.h"

#include <tudat/simulation/estimation.h>

#include "tudat/astro/gravitation/librationPoint.h"
#include "tudat/astro/gravitation/unitConversionsCircularRestrictedThreeBodyProblem.h"
#include "tudat/math/interpolators/createInterpolator.h"
#include "tudat/astro/basic_astro/celestialBodyConstants.h"

namespace tudat
{

namespace propagators
{

using namespace tudat::simulation_setup;

//! Function to create an integrator to propagate the dynamics (in normalized units) in CR3BP
std::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::Vector6d > > createCR3BPIntegrator(
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const double massParameter,
        const Eigen::Vector6d& initialState,
        const double initialTime )
{
    // Create state derivative model
    std::shared_ptr< StateDerivativeCircularRestrictedThreeBodyProblem > stateDerivativeModel =
            std::make_shared< StateDerivativeCircularRestrictedThreeBodyProblem  >( massParameter );
    std::function< Eigen::Vector6d( const double, const Eigen::Vector6d& ) > stateDerivativeFunction =
            std::bind( &StateDerivativeCircularRestrictedThreeBodyProblem::computeStateDerivative, stateDerivativeModel,
                       std::placeholders::_1, std::placeholders::_2 );

    // Create integrator object
    return numerical_integrators::createIntegrator< double, Eigen::Vector6d >(
                stateDerivativeFunction, initialState, initialTime, integratorSettings );
}

//! Function to propagate the dynamics (in normalized units) in CR3BP
std::map< double, Eigen::Vector6d > performCR3BPIntegration(
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const double massParameter,
        const Eigen::Vector6d& initialState,
        const double initialTime,
        const double finalTime,
        const bool propagateToExactFinalTime )
{
    // Create integrator object
    std::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::Vector6d > > integrator =
            createCR3BPIntegrator( integratorSettings, massParameter, initialState, initialTime );

    // Initialize return data map
    std::map< double, Eigen::Vector6d > stateHistory;

    // Store initial state and time
    double currentTime = initialTime;
    double secondToLastTime = initialTime;

    Eigen::Vector6d currentState = initialState;
    stateHistory[ currentTime ] = currentState;

    // Integrate to final time
    double timeStep = integratorSettings->initialTimeStep_;
    while( currentTime <= finalTime )
    {
        secondToLastTime = currentTime;
        currentState = integrator->performIntegrationStep( timeStep );
        currentTime = integrator->getCurrentIndependentVariable( );
        timeStep = integrator->getNextStepSize( );
        stateHistory[ currentTime ] = currentState;
    }

    if( propagateToExactFinalTime )
    {
        // Determine final time step and propagate
        double finalTimeStep = finalTime - secondToLastTime;

        stateHistory.erase( currentTime );
        integrator->rollbackToPreviousState( );

        currentState = integrator->performIntegrationStep( finalTimeStep );
        currentTime = integrator->getCurrentIndependentVariable( );
        stateHistory[ currentTime ] = currentState;
    }
    return stateHistory;

}


}

}
