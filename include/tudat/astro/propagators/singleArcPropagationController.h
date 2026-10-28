/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SINGLEARCPROPAGATIONCONTROLLER_H
#define TUDAT_SINGLEARCPROPAGATIONCONTROLLER_H

#include <chrono>
#include <map>

#include <Eigen/Core>

#include "tudat/basics/tudatTypeTraits.h"
#include "tudat/astro/propagators/dynamicsStateDerivativeModel.h"
#include "tudat/astro/propagators/integrateEquations.h"
#include "tudat/math/integrators/createNumericalIntegrator.h"
#include "tudat/simulation/propagation_setup/propagationProcessingSettings.h"
#include "tudat/simulation/propagation_setup/propagationTermination.h"

namespace tudat
{

namespace propagators
{

//! Class that links the state derivative model of a single-arc propagation to the numerical integrator.
/*!
 *  Class that links the state derivative model of a single-arc propagation (dynamics and, if requested, variational equations)
 *  to the numerical integrator. It holds the integrator and all objects needed during the propagation, and stores the
 *  numerical results as the propagation proceeds. Typical use is to call initialize, propagate and writeResults, in that order.
 *  \tparam StateType Type of the propagated state (column vector for dynamics only, matrix if variational equations are included)
 *  \tparam TimeType Type of the independent variable
 */
template< typename StateType = Eigen::MatrixXd, typename TimeType = double >
class SingleArcPropagationController
{
public:
    typedef typename StateType::Scalar StateScalarType;

    typedef typename scalar_type< TimeType >::value_type TimeStepType;

    //! Constructor
    /*!
     *  Constructor
     *  \param stateDerivativeModel Object that updates the environment and computes the state derivative
     *  \param integratorSettings Settings for the numerical integrator
     *  \param terminationCondition Object to determine when/how the propagation is to be stopped
     *  \param dependentVariableFunction Function returning dependent variables (empty if none are to be saved)
     *  \param processingSettings Settings for saving and printing of the propagation results
     */
    SingleArcPropagationController( const std::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > stateDerivativeModel,
                                    const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
                                    const std::shared_ptr< PropagationTerminationCondition > terminationCondition,
                                    const std::function< Eigen::VectorXd( ) > dependentVariableFunction,
                                    const std::shared_ptr< SingleArcPropagatorProcessingSettings > processingSettings ):
        stateDerivativeModel_( stateDerivativeModel ), integratorSettings_( integratorSettings ),
        terminationCondition_( terminationCondition ), dependentVariableFunction_( dependentVariableFunction ),
        processingSettings_( processingSettings )
    {
        statePostProcessingFunction_ =
                PostProcessingFunctionProvider< StateScalarType, TimeType, StateType::ColsAtCompileTime >::getPostProcessingFunction(
                        stateDerivativeModel );
        stateDerivativeFunction_ = [ = ]( const TimeType time, const StateType& state ) -> StateType {
            return stateDerivativeModel->computeStateDerivative( time, state );
        };
    }

    //! Function to create the numerical integrator, and save (and print, if requested) the initial conditions
    /*!
     *  Function to create the numerical integrator, and save (and print, if requested) the initial conditions
     *  \param initialState Initial state, in propagator-specific form
     *  \param initialTime Initial time
     */
    void initialize( const StateType& initialState, const TimeType initialTime )
    {
        integrator_ = numerical_integrators::createIntegrator< TimeType, StateType, TimeStepType >(
                stateDerivativeFunction_, initialState, initialTime, integratorSettings_ );
        if( integratorSettings_->assessTerminationOnMinorSteps_ )
        {
            integrator_->setPropagationTerminationFunction( std::bind( &PropagationTerminationCondition::checkStopCondition,
                                                                       terminationCondition_,
                                                                       std::placeholders::_1,
                                                                       std::placeholders::_2,
                                                                       std::placeholders::_3 ) );
        }

        initialClockTime_ = std::chrono::steady_clock::now( );
        terminationReason_ = std::make_shared< PropagationTerminationDetails >( unknown_propagation_termination_reason );

        // Add results at initial state
        currentTime_ = integrator_->getCurrentIndependentVariable( );
        initialTime_ = currentTime_;
        currentState_ = integrator_->getCurrentState( );

        solutionHistory_.clear( );
        solutionHistory_[ currentTime_ ] = currentState_;
        dependentVariableHistory_.clear( );
        if( dependentVariableFunction_ != nullptr )
        {
            updateEnvironment( );
            dependentVariableHistory_[ currentTime_ ] = dependentVariableFunction_( );
        }

        cumulativeComputationTimeHistory_.clear( );
        saveComputationTimeAndSignalEndOfMajorStep( );

        timeStep_ = integrator_->getNextStepSize( );
        stepsSinceLastSave_ = 1;
        timeOfLastSave_ = currentTime_;
        stepsSinceLastPrint_ = 1;
        timeOfLastPrint_ = TUDAT_NAN;

        if( processingSettings_->getPrintSettings( )->getPrintInitialAndFinalConditions( ) )
        {
            std::cout << "PRINTING INITIAL CONDITIONS" << std::endl;
            std::cout << "   Initial epoch: " << currentTime_ << std::endl;
            printState( "Initial state" );
            timeOfLastPrint_ = currentTime_;
            stepsSinceLastPrint_ = 0;
        }
    }

    //! Function to propagate the dynamics until the termination condition is reached
    void propagate( )
    {
        if( integrator_ == nullptr )
        {
            throw std::runtime_error( "Error when propagating dynamics, propagation controller is not initialized." );
        }

        propagateUntilCondition( terminationCondition_ );

        if( processingSettings_->getPrintSettings( )->getPrintInitialAndFinalConditions( ) )
        {
            std::cout << "PRINTING FINAL CONDITIONS" << std::endl;
            std::cout << "   Clock time since propagation start: " << currentCpuTime_ << std::endl;
            std::cout << "   Final epoch: " << currentTime_ << std::endl;
            printState( "Final state" );
        }
    }

    //! Function to set the numerical results of the propagation in a simulation results object
    template< typename SimulationResults >
    void writeResults( const std::shared_ptr< SimulationResults > simulationResults )
    {
        simulationResults->reset( solutionHistory_,
                                  dependentVariableHistory_,
                                  cumulativeComputationTimeHistory_,
                                  std::map< TimeType, unsigned int >( ),
                                  terminationReason_ );
    }

private:
    //! Function to propagate the dynamics until a given condition is met.
    /*!
     *  Function to propagate the dynamics until a given condition is met. Currently only called with the termination
     *  condition, but separated from it to allow the propagation to be paused at other conditions in the future.
     *  \param condition Condition at which propagation is to be stopped
     */
    void propagateUntilCondition( const std::shared_ptr< PropagationTerminationCondition > condition )
    {
        bool breakPropagation = false;
        do
        {
            try
            {
                bool environmentIsUpdated = false;
                if( currentState_.allFinite( ) && !currentState_.hasNaN( ) )
                {
                    printCurrentStateIfRequired( );
                    if( performIntegrationStep( ) )
                    {
                        // Termination condition reached during integrator sub-steps; last step not computed correctly, so
                        // return immediately without saving the current state.
                        terminationReason_ = std::make_shared< PropagationTerminationDetails >( termination_condition_reached, 0 );
                        break;
                    }
                    environmentIsUpdated = saveCurrentStepIfRequired( );
                }
                else
                {
                    std::cerr << "Error, propagation terminated at t=" + std::to_string( static_cast< double >( currentTime_ ) ) +
                                    ", found NaN/Inf entry, returning propagation data up to current time. Current state "
                                    "(transposed) is: "
                              << currentState_.transpose( ) << std::endl;
                    breakPropagation = true;
                    terminationReason_ = std::make_shared< PropagationTerminationDetails >( nan_or_inf_detected_in_state );
                }

                if( !environmentIsUpdated && condition->requiresEnvironmentUpdate( ) )
                {
                    updateEnvironment( );
                }

                if( checkStopCondition( condition ) )
                {
                    breakPropagation = true;
                }

                saveComputationTimeAndSignalEndOfMajorStep( );
            }
            catch( const std::exception& caughtException )
            {
                std::cerr << caughtException.what( ) << std::endl;
                std::cerr << "Error, propagation terminated at t=" + std::to_string( static_cast< double >( currentTime_ ) ) +
                                ", returning propagation data up to current time."
                          << std::endl;
                breakPropagation = true;
                terminationReason_ = std::make_shared< PropagationTerminationDetails >( runtime_error_caught_in_propagation );
            }
        } while( !breakPropagation );
    }

    //! Function to perform a single integration step, and post-process the new state.
    /*!
     *  Function to perform a single integration step, and post-process the new state.
     *  \return True if the termination condition was reached during the integrator sub-steps (in which case the new
     *  state is not valid)
     */
    bool performIntegrationStep( )
    {
        currentState_ = integrator_->performIntegrationStep( timeStep_ );
        statePostProcessingFunction_( currentState_ );
        integrator_->modifyCurrentState( currentState_, true );

        // Only returns true if evaluation of the termination condition during integrator sub-steps is enabled.
        if( integrator_->getPropagationTerminationConditionReached( ) )
        {
            return true;
        }

        currentTime_ = integrator_->getCurrentIndependentVariable( );
        timeStep_ = integrator_->getNextStepSize( );
        return false;
    }

    //! Function to update the environment and state derivative model to the current time and state.
    void updateEnvironment( )
    {
        integrator_->getStateDerivativeFunction( )( currentTime_, currentState_ );
    }

    //! Function to save the current state and dependent variables, if required by the processing settings.
    /*!
     *  Function to save the current state and dependent variables, if required by the processing settings.
     *  \return True if the environment was updated to the current state (to compute the dependent variables)
     */
    bool saveCurrentStepIfRequired( )
    {
        bool environmentIsUpdated = false;
        if( processingSettings_->saveCurrentStep( stepsSinceLastSave_,
                                                  std::fabs( static_cast< double >( currentTime_ ) - timeOfLastSave_ ) ) )
        {
            solutionHistory_[ currentTime_ ] = currentState_;
            if( dependentVariableFunction_ != nullptr )
            {
                updateEnvironment( );
                environmentIsUpdated = true;
                dependentVariableHistory_[ currentTime_ ] = dependentVariableFunction_( );
            }
            timeOfLastSave_ = currentTime_;
            stepsSinceLastSave_ = 0;
        }
        stepsSinceLastSave_++;
        return environmentIsUpdated;
    }

    //! Function to print the current state, if required by the print settings.
    void printCurrentStateIfRequired( )
    {
        if( processingSettings_->getPrintSettings( )->printCurrentStep(
                    stepsSinceLastPrint_, std::fabs( static_cast< double >( currentTime_ ) - timeOfLastPrint_ ) ) )
        {
            stepsSinceLastPrint_ = 0;
            timeOfLastPrint_ = currentTime_;
            std::cout << "PRINTING STATE DURING PROPAGATION" << std::endl;
            std::cout << "   Clock time since propagation start: " << currentCpuTime_ << std::endl;
            std::cout << "   Time since initial epoch: " << currentTime_ - initialTime_ << std::endl;
            printState( "Current state" );
        }
        stepsSinceLastPrint_++;
    }

    //! Function to check whether the propagation is to be stopped and, if so, propagate to the exact condition if required.
    /*!
     *  Function to check whether the propagation is to be stopped and, if so, propagate to the exact condition if required
     *  and set the termination details.
     *  \param condition Condition at which propagation is to be stopped
     *  \return True if propagation is to be stopped
     */
    bool checkStopCondition( const std::shared_ptr< PropagationTerminationCondition > condition )
    {
        if( !condition->checkStopCondition(
                    static_cast< double >( currentTime_ ), currentCpuTime_, currentState_.template cast< double >( ) ) )
        {
            return false;
        }

        if( condition->iterateToExactTermination( ) )
        {
            currentTime_ = propagateToExactTerminationCondition( integrator_,
                                                                 condition,
                                                                 timeStep_,
                                                                 dependentVariableFunction_,
                                                                 solutionHistory_,
                                                                 dependentVariableHistory_,
                                                                 currentCpuTime_ );
        }

        if( condition->getTerminationType( ) != hybrid_stopping_condition )
        {
            terminationReason_ = std::make_shared< PropagationTerminationDetails >( termination_condition_reached,
                                                                                    condition->getcheckTerminationToExactCondition( ) );
        }
        else
        {
            std::shared_ptr< HybridPropagationTerminationCondition > hybridCondition =
                    std::dynamic_pointer_cast< HybridPropagationTerminationCondition >( condition );
            if( hybridCondition == nullptr )
            {
                throw std::runtime_error( "Error when saving termination reason, type is hybrid, but class is not." );
            }
            terminationReason_ = std::make_shared< PropagationTerminationDetailsFromHybridCondition >(
                    condition->iterateToExactTermination( ), hybridCondition );
        }
        return true;
    }

    //! Function to save the cumulative computation time at the current time, and signal end of step to the state derivative model.
    void saveComputationTimeAndSignalEndOfMajorStep( )
    {
        stateDerivativeModel_->signalEndOfMajorStep( currentTime_ );
        currentCpuTime_ =
                std::chrono::duration_cast< std::chrono::nanoseconds >( std::chrono::steady_clock::now( ) - initialClockTime_ ).count( ) *
                1.0e-9;
        cumulativeComputationTimeHistory_[ currentTime_ ] = currentCpuTime_;
    }

    //! Function to print the current state to the terminal, with a given label.
    void printState( const std::string& label )
    {
        if( currentState_.cols( ) == 1 )
        {
            std::cout << "   " << label << " (transpose): " << std::endl << currentState_.transpose( ) << std::endl << std::endl;
        }
        else
        {
            std::cout << "   " << label << ": " << std::endl << currentState_ << std::endl << std::endl;
        }
    }

    //! Object that updates the environment and computes the state derivative
    std::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > stateDerivativeModel_;

    //! Function computing the state derivative from the current time and state (calls stateDerivativeModel_).
    std::function< StateType( const TimeType, const StateType& ) > stateDerivativeFunction_;

    //! Settings for the numerical integrator
    std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings_;

    //! Object to determine when/how the propagation is to be stopped
    std::shared_ptr< PropagationTerminationCondition > terminationCondition_;

    //! Function returning dependent variables (empty if none are to be saved)
    std::function< Eigen::VectorXd( ) > dependentVariableFunction_;

    //! Function to post-process state after each integration step (obtained from stateDerivativeModel_)
    std::function< void( StateType& ) > statePostProcessingFunction_;

    //! Settings for saving and printing of the propagation results
    std::shared_ptr< SingleArcPropagatorProcessingSettings > processingSettings_;

    //! Numerical integrator, created by initialize
    std::shared_ptr< numerical_integrators::NumericalIntegrator< TimeType, StateType, StateType, TimeStepType > > integrator_;

    //! Saved history of numerically integrated states
    std::map< TimeType, StateType > solutionHistory_;

    //! Saved history of dependent variables
    std::map< TimeType, Eigen::VectorXd > dependentVariableHistory_;

    //! History of cumulative computation time
    std::map< TimeType, double > cumulativeComputationTimeHistory_;

    //! Reason for which the propagation was terminated
    std::shared_ptr< PropagationTerminationDetails > terminationReason_;

    //! Clock time at start of propagation
    std::chrono::steady_clock::time_point initialClockTime_;

    //! Initial time of propagation
    TimeType initialTime_;

    //! Current time of propagation
    TimeType currentTime_;

    //! Current propagated state
    StateType currentState_;

    //! Step size to be used for next integration step
    TimeStepType timeStep_;

    //! Computation time since start of propagation
    double currentCpuTime_ = 0.0;

    //! Number of steps since state was last saved
    int stepsSinceLastSave_;

    //! Time at which state was last saved
    double timeOfLastSave_;

    //! Number of steps since state was last printed
    int stepsSinceLastPrint_;

    //! Time at which state was last printed
    double timeOfLastPrint_;
};

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_SINGLEARCPROPAGATIONCONTROLLER_H
