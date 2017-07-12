/*    Copyright (c) 2010-2017, Delft University of Technology
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

boost::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::Vector6d > > createCR3BPIntegrator(
        const boost::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const double massParameter,
        const Eigen::Vector6d& initialState )
{
    boost::shared_ptr< StateDerivativeCircularRestrictedThreeBodyProblem > stateDerivativeModel =
            boost::make_shared< StateDerivativeCircularRestrictedThreeBodyProblem  >( massParameter );

    boost::function< Eigen::Vector6d( const double, const Eigen::Vector6d& ) > stateDerivativeFunction =
            boost::bind( &StateDerivativeCircularRestrictedThreeBodyProblem::computeStateDerivative, stateDerivativeModel,
                         _1, _2 );
    return numerical_integrators::createIntegrator< double, Eigen::Vector6d >(
                stateDerivativeFunction, initialState, integratorSettings );
}

std::map< double, Eigen::Vector6d > performCR3BPIntegration(
        const boost::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const double massParameter,
        const Eigen::Vector6d& initialState,
        const double finalTime )
{
    boost::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::Vector6d > > integrator =
            createCR3BPIntegrator( integratorSettings, massParameter, initialState );
    std::map< double, Eigen::Vector6d > stateHistory;

    double currentTime = integratorSettings->initialTime_;
    Eigen::Vector6d currentState = initialState;
    stateHistory[ currentTime ] = currentState;

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

}

}
