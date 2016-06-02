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

namespace tudat
{

namespace propagators
{

template< typename OutputType, typename InputType >
OutputType evaluateReferenceFunction(
        const boost::function< OutputType( const InputType&, const InputType& ) > functionToEvaluate,
        const boost::function< InputType( ) > firstInput,
        const boost::function< InputType( ) > secondInput )
{
    return functionToEvaluate( firstInput( ), secondInput( ) );
}

template< typename OutputType, typename InputType >
OutputType evaluateFunction(
        const boost::function< OutputType( const InputType, const InputType ) > functionToEvaluate,
        const boost::function< InputType( ) > firstInput,
        const boost::function< InputType( ) > secondInput )
{
    return functionToEvaluate( firstInput( ), secondInput( ) );
}

boost::function< double( ) > getDependentVariableFunction(
        const PropagationDependentVariables dependentVariable,
        const std::string& bodyWithProperty,
        const std::string& secondaryBody,
        const simulation_setup::NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationModelList = basic_astrodynamics::AccelerationMap( ) );

class PropagationPropagationTerminationCondition
{
public:
    PropagationPropagationTerminationCondition( ){ }

    virtual ~PropagationPropagationTerminationCondition( ){ }

    virtual bool checkStopCondition( const double time ) = 0;
};

class FixedTimePropagationTerminationCondition: public PropagationPropagationTerminationCondition
{
public:
    FixedTimePropagationTerminationCondition(
            const double stopTime,
            const bool propagationDirectionIsPositive ):
        stopTime_( stopTime ),
        propagationDirectionIsPositive_( propagationDirectionIsPositive ){ }


    bool checkStopCondition( const double time );

private:
    double stopTime_;

    bool propagationDirectionIsPositive_;
};

class SingleVariableLimitPropagationTerminationCondition: public PropagationPropagationTerminationCondition
{
public:
    SingleVariableLimitPropagationTerminationCondition(
            const std::pair< PropagationDependentVariables, std::string > variableType,
            const boost::function< double( ) > variableRetrievalFuntion,
            const double limitingValue,
            const bool useAsLowerBound ):
    variableType_( variableType ), variableRetrievalFuntion_( variableRetrievalFuntion ),
    limitingValue_( limitingValue ), useAsLowerBound_( useAsLowerBound ){ }

    virtual ~SingleVariableLimitPropagationTerminationCondition( ){ }

    bool checkStopCondition( const double time );

private:
    std::pair< PropagationDependentVariables, std::string > variableType_;

    boost::function< double( ) > variableRetrievalFuntion_;

    double limitingValue_;

    bool useAsLowerBound_;
};

class HybridPropagationTerminationCondition: public PropagationPropagationTerminationCondition
{
public:
    HybridPropagationTerminationCondition(
            const std::vector< boost::shared_ptr< PropagationPropagationTerminationCondition > > propagationTerminationCondition,
            const bool fulFillSingleCondition = 0 ):
        propagationTerminationCondition_( propagationTerminationCondition ), fulFillSingleCondition_( fulFillSingleCondition ){ }

    bool checkStopCondition( const double time );

private:

    std::vector< boost::shared_ptr< PropagationPropagationTerminationCondition > > propagationTerminationCondition_;

    bool fulFillSingleCondition_;
};


boost::shared_ptr< PropagationPropagationTerminationCondition > createPropagationPropagationTerminationConditions(
        const boost::shared_ptr< PropagationTerminationSettings > terminationSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const double initialTimeStep );

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
std::map< TimeType, StateType > integrateEquations(
        boost::function< StateType( const TimeType, const StateType&) > stateDerivativeFunction,
        const StateType initialState,
        const boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
        const boost::function< bool( const double ) > stopPropagationFunction,
        const TimeType printInterval = TUDAT_NAN )
{
    using namespace tudat::numerical_integrators;


    // Create numerical integrator.
    boost::shared_ptr< NumericalIntegrator< TimeType, StateType, StateType > > integrator =
            createIntegrator< TimeType, StateType >( stateDerivativeFunction, initialState, integratorSettings );

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
    TimeType previousTime = currentTime;

    // Perform first integration step.
    newState = integrator->performIntegrationStep( timeStep );

    currentTime = integrator->getCurrentIndependentVariable( );

    timeStep = timeStepSign * integrator->getNextStepSize( );
    solutionHistory[ currentTime ] = newState;

    int printIndex = 0;
    int printFrequency = integratorSettings->printFrequency_;
    // Perform numerical integration steps until end time reached.
    do
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
                std::cout<<"Current time and state in integration: "<<std::setprecision( 10 )<<
                           timeStep<<" "<<currentTime<<" "<<newState.transpose( )<<std::endl;
            }
        }
    }
    while( !stopPropagationFunction( static_cast< double >( currentTime ) ) );

    return solutionHistory;
}

}

}
#endif // TUDAT_INTEGRATEEQUATIONS_H
