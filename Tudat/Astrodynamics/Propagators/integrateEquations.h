/*    Copyright (c) 2010-2017, Delft University of Technology
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
#include <boost/lambda/lambda.hpp>

#include <map>

#include "Tudat/Mathematics/NumericalIntegrators/numericalIntegrator.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Basics/timeType.h"
#include "Tudat/Astrodynamics/Propagators/singleStateTypeDerivative.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationTermination.h"

namespace tudat
{

namespace propagators
{

//! Function to numerically integrate a given first order differential equation
/*!
 *  Function to numerically integrate a given first order differential equation, with the state derivative a function of
 *  a single independent variable and the current state
 *  \param integrator Numerical integrator used for propagation
 *  \param initialTimeStep Time step to use for first step of numerical integration
 *  \param stopPropagationFunction Function determining whether the propagation is to be stopped at the current time.
 *  \param solutionHistory History of dependent variables that are to be saved given as map
 *  (time as key; returned by reference)
 *  \param dependentVariableHistory History of dependent variables that are to be saved given as map
 *  (time as key; returned by reference)
 *  \param propagationTerminationReason Event that triggered the termination of the propagation (returned by reference)
 *  \param dependentVariableFunction Function returning dependent variables (obtained from environment and state
 *  derivative model).
 *  \param saveFrequency Frequency at which to save the numerical integrated states (in units of i.e. per n integration time
 *  steps, with n = saveFrequency).
 *  \param printInterval Frequency with which to print progress to console (nan = never).
 */
template< typename StateType = Eigen::MatrixXd, typename TimeType = double, typename TimeStepType = TimeType  >
void integrateEquationsFromIntegrator(
        const boost::shared_ptr< numerical_integrators::NumericalIntegrator< TimeType, StateType, StateType, TimeStepType > > integrator,
        const TimeStepType initialTimeStep,
        const boost::function< bool( const double ) > stopPropagationFunction,
        std::map< TimeType, StateType >& solutionHistory,
        std::map< TimeType, Eigen::VectorXd >& dependentVariableHistory,
        PropagationTerminationReason &propagationTerminationReason,
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


    if( !dependentVariableFunction.empty( ) )
    {
        integrator->getStateDerivativeFunction( )( currentTime, newState );
        dependentVariableHistory[ currentTime ] = dependentVariableFunction( );
    }

    // Set initial time step and total integration time.
    TimeStepType timeStep = initialTimeStep;
    TimeType previousTime = currentTime;

    int saveIndex = 0;

    propagationTerminationReason = termination_condition_reached;
    bool breakPropagation = 0;
    // Perform numerical integration steps until end time reached.
    do
    {
        try
        {
            previousTime = currentTime;

            // Perform integration step.
            newState = integrator->performIntegrationStep( timeStep );
            currentTime = integrator->getCurrentIndependentVariable( );
            timeStep = integrator->getNextStepSize( );

            // Check if the termination condition was reached during evaluation of integration substeps.
            // If evaluation of the termination during integration substeps is disabled (default behaviour for all
            // propagators execpt `dsst`) this function returns always `false`.
            breakPropagation = integrator->getPropagationShouldTerminate();

            if ( !breakPropagation )
            {
                // Save integration result in map
                saveIndex++;
                saveIndex = saveIndex % saveFrequency;
                if( saveIndex == 0 )
                {
                    solutionHistory[ currentTime ] = newState;

                    if( !dependentVariableFunction.empty( ) )
                    {
                        integrator->getStateDerivativeFunction( )( currentTime, newState );
                        dependentVariableHistory[ currentTime ] = dependentVariableFunction( );
                    }
                }


                // Print solutions
                if( printInterval == printInterval && !breakPropagation )
                {
                    if( ( static_cast<int>( std::fabs( static_cast< double >( currentTime - initialTime ) ) ) %
                          static_cast< int >( printInterval ) ) <
                            ( static_cast< int >( std::fabs( static_cast< double >( previousTime - initialTime ) ) ) %
                              static_cast<int>( printInterval ) )  )
                    {
                        std::cout<<"Current time and state in integration: "<<std::setprecision( 10 )<<
                                   timeStep<<" "<<currentTime<<" "<<newState.transpose( )<<std::endl;
                    }
                }
            }
        }
        catch( const std::exception &caughtException )
        {
            std::cerr << caughtException.what( )<<std::endl;
            std::cerr<<"Error, propagation terminated at t=" + boost::lexical_cast< std::string >( currentTime ) +
                       ", returning propagation data up to current time"<<std::endl;
            breakPropagation = 1;
            propagationTerminationReason = runtime_error_caught;
        }
    }
    while( !stopPropagationFunction( static_cast< double >( currentTime ) ) && !breakPropagation );
}


//! Interface class for integrating some state derivative function.
/*!
 *  Interface class for integrating some state derivative function.. This class is used instead of a single templated free
 *  function to allow ObservationModel the integrator etc. to adapt its time step variable to long double if the Time
 *  object is used as TimeType. This class has template specializations for double/Time TimeType, and contains a single
 *  integrateEquations function that performs the required operation.
 */
template< typename StateType = Eigen::MatrixXd, typename TimeType = double >
class EquationIntegrationInterface
{
public:

    //! Function to numerically integrate a given first order differential equation
    /*!
     *  Function to numerically integrate a given first order differential equation, with the state derivative a function of
     *  a single independent variable and the current state
     *  \param stateDerivativeFunction Function returning the state derivative from current time and state.
     *  \param solutionHistory History of numerical states given as map (time as key; returned by reference)
     *  \param initialState Initial state
     *  \param integratorSettings Settings for numerical integrator.
     *  \param stopPropagationFunction Function determining whether the propagation is to be stopped at the current time.
     *  \param dependentVariableHistory History of dependent variables that are to be saved given as map
     *  (time as key; returned by reference)
     *  \param propagationTerminationReason Event that triggered the termination of the propagation (returned by
     *  reference)
     *  \param dependentVariableFunction Function returning dependent variables (obtained from environment and state
     *  derivative model).
     *  \param printInterval Frequency with which to print progress to console (nan = never).
     *  \param assessPropagationTerminationConditionDuringIntegrationSubsteps Whether the propagation termination
     *  conditions should be evaluated during the intermediate state updates performed by the integrator to compute
     *  the quantities necessary to integrate the state to a new epoch (default value is false).
     */
    static void integrateEquations(
            boost::function< StateType( const TimeType, const StateType& ) > stateDerivativeFunction,
            std::map< TimeType, StateType >& solutionHistory,
            const StateType initialState,
            const boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const boost::function< bool( const double ) > stopPropagationFunction,
            std::map< TimeType, Eigen::VectorXd >& dependentVariableHistory,
            PropagationTerminationReason &propagationTerminationReason,
            const boost::function< Eigen::VectorXd( ) > dependentVariableFunction =
            boost::function< Eigen::VectorXd( ) >( ),
            const TimeType printInterval = TUDAT_NAN,
            const bool assessPropagationTerminationConditionDuringIntegrationSubsteps = false );
};

//! Interface class for integrating some state derivative function.
template< typename StateType >
class EquationIntegrationInterface< StateType, double >
{
public:

    //! Function to numerically integrate a given first order differential equation
    /*!
     *  Function to numerically integrate a given first order differential equation, with the state derivative a function of
     *  a single independent variable and the current state
     *  \param stateDerivativeFunction Function returning the state derivative from current time and state.
     *  \param solutionHistory History of numerical states given as map (time as key; returned by reference)
     *  \param initialState Initial state
     *  \param integratorSettings Settings for numerical integrator.
     *  \param stopPropagationFunction Function determining whether the propagation is to be stopped at the current time.
     *  \param dependentVariableHistory History of dependent variables that are to be saved given as map
     *  (time as key; returned by reference)
     *  \param propagationTerminationReason Event that triggered the termination of the propagation (returned by
     *  reference)
     *  \param dependentVariableFunction Function returning dependent variables (obtained from environment and state
     *  derivative model).
     *  \param printInterval Frequency with which to print progress to console (nan = never).
     *  \param assessPropagationTerminationConditionDuringIntegrationSubsteps Whether the propagation termination
     *  conditions should be evaluated during the intermediate state updates performed by the integrator to compute
     *  the quantities necessary to integrate the state to a new epoch (default value is false).
     */
    static void integrateEquations(
            boost::function< StateType( const double, const StateType& ) > stateDerivativeFunction,
            std::map< double, StateType >& solutionHistory,
            const StateType initialState,
            const boost::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
            const boost::function< bool( const double ) > stopPropagationFunction,
            std::map< double, Eigen::VectorXd >& dependentVariableHistory,
            PropagationTerminationReason &propagationTerminationReason,
            const boost::function< Eigen::VectorXd( ) > dependentVariableFunction =
            boost::function< Eigen::VectorXd( ) >( ),
            const double printInterval = TUDAT_NAN,
            const bool assessPropagationTerminationConditionDuringIntegrationSubsteps = false )
    {
        // Create numerical integrator.
        boost::shared_ptr< numerical_integrators::NumericalIntegrator< double, StateType, StateType > > integrator =
                numerical_integrators::createIntegrator< double, StateType >(
                    stateDerivativeFunction, initialState, integratorSettings );

        if ( assessPropagationTerminationConditionDuringIntegrationSubsteps )
        {
            integrator->setPropagationTerminationFunction( stopPropagationFunction );
        }

        integrateEquationsFromIntegrator< StateType, double >(
                    integrator, integratorSettings->initialTimeStep_, stopPropagationFunction, solutionHistory,
                    dependentVariableHistory, propagationTerminationReason,
                    dependentVariableFunction,
                    integratorSettings->saveFrequency_, printInterval );
    }
};

//! Interface class for integrating some state derivative function.
template< typename StateType >
class EquationIntegrationInterface< StateType, Time >
{
public:

    //! Function to numerically integrate a given first order differential equation
    /*!
     *  Function to numerically integrate a given first order differential equation, with the state derivative a function of
     *  a single independent variable and the current state
     *  \param stateDerivativeFunction Function returning the state derivative from current time and state.
     *  \param solutionHistory History of numerical states given as map (time as key; returned by reference)
     *  \param initialState Initial state
     *  \param integratorSettings Settings for numerical integrator.
     *  \param stopPropagationFunction Function determining whether the propagation is to be stopped at the current time.
     *  \param dependentVariableHistory History of dependent variables that are to be saved given as map
     *  (time as key; returned by reference)
     *  \param propagationTerminationReason Event that triggered the termination of the propagation (returned by
     *  reference)
     *  \param dependentVariableFunction Function returning dependent variables (obtained from environment and state
     *  derivative model).
     *  \param printInterval Frequency with which to print progress to console (nan = never).
     *  \param assessPropagationTerminationConditionDuringIntegrationSubsteps Whether the propagation termination
     *  conditions should be evaluated during the intermediate state updates performed by the integrator to compute
     *  the quantities necessary to integrate the state to a new epoch (default value is false).
     */
    static void integrateEquations(
            boost::function< StateType( const Time, const StateType& ) > stateDerivativeFunction,
            std::map< Time, StateType >& solutionHistory,
            const StateType initialState,
            const boost::shared_ptr< numerical_integrators::IntegratorSettings< Time > > integratorSettings,
            const boost::function< bool( const double ) > stopPropagationFunction,
            std::map< Time, Eigen::VectorXd >& dependentVariableHistory,
            PropagationTerminationReason &propagationTerminationReason,
            const boost::function< Eigen::VectorXd( ) > dependentVariableFunction =
            boost::function< Eigen::VectorXd( ) >( ),
            const Time printInterval = TUDAT_NAN,
            const bool assessPropagationTerminationConditionDuringIntegrationSubsteps = false )
    {
        // Create numerical integrator.
        boost::shared_ptr< numerical_integrators::NumericalIntegrator< Time, StateType, StateType, long double > > integrator =
                numerical_integrators::createIntegrator< Time, StateType, long double  >(
                    stateDerivativeFunction, initialState, integratorSettings );

        if ( assessPropagationTerminationConditionDuringIntegrationSubsteps )
        {
            integrator->setPropagationTerminationFunction( stopPropagationFunction );
        }

        integrateEquationsFromIntegrator< StateType, Time, long double >(
                    integrator, integratorSettings->initialTimeStep_, stopPropagationFunction, solutionHistory,
                    dependentVariableHistory, propagationTerminationReason,
                    dependentVariableFunction,
                    integratorSettings->saveFrequency_, printInterval );
    }
};

} // namespace propagators

} // namespace tudat

#endif // TUDAT_INTEGRATEEQUATIONS_H
