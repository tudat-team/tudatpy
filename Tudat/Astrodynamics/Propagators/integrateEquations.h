/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_INTEGRATEEQUATIONS_H
#define TUDAT_INTEGRATEEQUATIONS_H

#include <Eigen/Core>

#include <map>

#include "Tudat/Mathematics/NumericalIntegrators/numericalIntegrator.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/Propagators/singleStateTypeDerivative.h"
#include "Tudat/Astrodynamics/Propagators/environmentUpdater.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"
#include "Tudat/Astrodynamics/Propagators/propagationTerminationConditions.h"

namespace tudat
{

namespace propagators
{


template< typename StateType = Eigen::MatrixXd, typename TimeType = double >
void integrateEquations(
        const boost::shared_ptr< numerical_integrators::NumericalIntegrator< TimeType, StateType, StateType > > integrator,
        const double initialTimeStep,
        const boost::function< bool( const double ) > stopPropagationFunction,
        std::map< TimeType, StateType >& solutionHistory,
        std::map< TimeType, Eigen::VectorXd >& dependentVariableHistory,
        const boost::function< Eigen::VectorXd( ) > dependentVariableFunction =
        boost::function< Eigen::VectorXd( ) >( ),
        const int saveFrequency = TUDAT_NAN,
        const TimeType printInterval = TUDAT_NAN )
{

    // Get Initial state and time.
    TimeType currentTime = integrator->getCurrentIndependentVariable( );
    TimeType initialTime = currentTime;
    StateType newState = integrator->getCurrentState( );

    // Initialization of numerical solutions for variational equations
    solutionHistory.clear( );
    solutionHistory[ currentTime ] = newState;
    dependentVariableHistory.clear( );

    // Check if numerical integration is forward or backwrd.
    TimeType timeStepSign = 1.0L;
//    if( initialTimeStep < 0.0 )
//    {
//        timeStepSign = -1.0L;
//    }

    // Set initial time step and total integration time.
    TimeType timeStep = initialTimeStep;
    TimeType previousTime = currentTime;

    // Perform first integration step.
    newState = integrator->performIntegrationStep( timeStep );

    currentTime = integrator->getCurrentIndependentVariable( );

    timeStep = timeStepSign * integrator->getNextStepSize( );
    solutionHistory[ currentTime ] = newState;
    if( !dependentVariableHistory.empty( ) )
    {
        dependentVariableHistory[ currentTime ] = dependentVariableFunction( );
    }

    int saveIndex = 0;

    // Perform numerical integration steps until end time reached.
    do
    {
        previousTime = currentTime;

        // Perform integration step.
        newState = integrator->performIntegrationStep( timeStep );
        currentTime = integrator->getCurrentIndependentVariable( );
        timeStep = timeStepSign * integrator->getNextStepSize( );

        // Save integration result in map
        saveIndex++;
        saveIndex = saveIndex % saveFrequency;
        if( saveIndex == 0 )
        {
            solutionHistory[ currentTime ] = newState;

            if( !dependentVariableHistory.empty( ) )
            {
                dependentVariableHistory[ currentTime ] = dependentVariableFunction( );
            }
        }

        // Print solutions
        if( printInterval == printInterval )
        {
            if( ( static_cast<int>( std::fabs( currentTime - initialTime ) ) %
                  static_cast< int >( printInterval ) ) <
                    ( static_cast< int >( std::fabs( previousTime - initialTime ) ) %
                      static_cast<int>( printInterval ) )  )
            {
                std::cout<<"Current time and state in integration: "<<std::setprecision( 10 )<<
                           timeStep<<" "<<currentTime<<" "<<newState.transpose( )<<std::endl;
            }
        }
    }
    while( !stopPropagationFunction( static_cast< double >( currentTime ) ) );
}

//! Function to numerically integrate a given first order differential equation
/*!
 *  Function to numerically integrate a given first order differential equation, with the state derivative a function of
 *  a single independent variable and the current state
 *  \param stateDerivativeFunction Function returning the state derivative from current time and state.
 *  \param initialState Initial state
 *  \param integratorSettings Settings for numerical integrator.
 *  \param printInterval Frequency with which to print progress to console (nan = never).
 *  \return History of numerical states (first of pair) and derivatives of states (second of pair) given as maps with time
 *  as key.
 */
template< typename StateType = Eigen::MatrixXd, typename TimeType = double >
void integrateEquations(
        boost::function< StateType( const TimeType, const StateType& ) > stateDerivativeFunction,
        std::map< TimeType, StateType >& solutionHistory,
        const StateType initialState,
        const boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
        const boost::function< bool( const double ) > stopPropagationFunction,
        const boost::function< Eigen::VectorXd( ) > dependentVariableFunction =
        boost::function< Eigen::VectorXd( ) >( ),
        const TimeType printInterval = TUDAT_NAN )
{
    std::map< TimeType, Eigen::VectorXd > dependentVariableHistory;

    // Create numerical integrator.
    boost::shared_ptr< numerical_integrators::NumericalIntegrator< TimeType, StateType, StateType > > integrator =
            numerical_integrators::createIntegrator< TimeType, StateType >(
                stateDerivativeFunction, initialState, integratorSettings );

    integrateEquations< StateType, TimeType >(
                integrator, integratorSettings->initialTimeStep_, stopPropagationFunction, solutionHistory,
                dependentVariableHistory,
                dependentVariableFunction,
                integratorSettings->saveFrequency_, printInterval );

}

}

}
#endif // TUDAT_INTEGRATEEQUATIONS_H
