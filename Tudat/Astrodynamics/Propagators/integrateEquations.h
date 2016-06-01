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

#include <iostream>
#include <map>

#include <Eigen/Core>

#include "Tudat/Mathematics/NumericalIntegrators/numericalIntegrator.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/Propagators/singleStateTypeDerivative.h"
#include "Tudat/Astrodynamics/Propagators/environmentUpdater.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"

namespace tudat
{

namespace propagators
{

//! Function to numerically integrate a given first order differential equation
/*!
 *  Function to numerically integrate a given first order differential equation, with the state
 *  derivative a function of a single independent variable and the current state
 *  \param stateDerivativeFunction Function returning the state derivative from current time and state.
 *  \param initialState Initial state
 *  \param integratorSettings Settings for numerical integrator.
 *  \param printInterval Frequency with which to print progress to console (nan = never).
 *  \return History of numerical states (first of pair) and derivatives of states (second of pair)
 *  given as maps with time as key.
 */
template< typename StateType = Eigen::MatrixXd, typename TimeType = double >
std::map< TimeType, StateType > integrateEquations(
        boost::function< StateType( const TimeType, const StateType&) > stateDerivativeFunction,
        const StateType initialState,
        boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
        const TimeType printInterval = TUDAT_NAN )
{
    using namespace tudat::numerical_integrators;


    // Create numerical integrator.
    boost::shared_ptr< NumericalIntegrator< TimeType, StateType, StateType > > integrator
          = createIntegrator< TimeType, StateType >( stateDerivativeFunction,
                                                     initialState, integratorSettings );

    // Get Initial state and time.
    TimeType currentTime = integratorSettings->initialTime_;
    StateType newState = initialState;

    // Initialization of numerical solutions for variational equations.
    std::map< TimeType, StateType > solutionHistory;
    solutionHistory[ currentTime ] = newState;

    // Check if numerical integration is forward or backwrd.
    TimeType timeStepSign = 1.0L;
    if( integratorSettings->initialTimeStep_ < 0.0 )
    {
        timeStepSign = -1.0L;
    }

    // Set initial time step and total integration time.
    TimeType timeStep = integratorSettings->initialTimeStep_;
    TimeType endTime = integratorSettings->endTime_;
    TimeType previousTime = currentTime;

    // Perform first integration step.
    newState = integrator->performIntegrationStep( timeStep );

    currentTime = integrator->getCurrentIndependentVariable( );

    timeStep = timeStepSign * integrator->getNextStepSize( );
    solutionHistory[ currentTime ] = newState;

    int printIndex = 0;
    int printFrequency = integratorSettings->printFrequency_;

    // Perform numerical integration steps until end time reached.
    while( timeStepSign * static_cast< TimeType >( currentTime )
           < timeStepSign * static_cast< TimeType >( endTime ) )
    {
        previousTime = currentTime;

        // Perform integration step.
        newState = integrator->performIntegrationStep( timeStep );
        currentTime = integrator->getCurrentIndependentVariable( );
        timeStep = timeStepSign * integrator->getNextStepSize( );

        // Save integration result in map
        printIndex++;
        printIndex = printIndex % printFrequency;
        if( printIndex == 0 )
        {
            solutionHistory[ currentTime ] = newState;
        }

        // Print solutions
        if( printInterval == printInterval )
        {
            if( ( static_cast<int>( std::fabs( currentTime - integratorSettings->initialTime_ ) ) %
                  static_cast< int >( printInterval ) ) <
                    ( static_cast< int >( std::fabs( previousTime - integratorSettings->initialTime_ ) ) %
                      static_cast<int>( printInterval ) )  )
            {
                std::cout << "Current time and state in integration: " << std::setprecision( 10 ) <<
                           timeStep << " " << currentTime << " " << newState.transpose( ) << std::endl;
            }
        }
    }

    return solutionHistory;
}

} // namespace propagators

} // namespace tudat
#endif // TUDAT_INTEGRATEEQUATIONS_H
