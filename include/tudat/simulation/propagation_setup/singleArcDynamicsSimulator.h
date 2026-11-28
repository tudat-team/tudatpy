/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SINGLEARCDYNAMICSSIMULATOR_H
#define TUDAT_SINGLEARCDYNAMICSSIMULATOR_H

#include "tudat/simulation/propagation_setup/dynamicsSimulatorBase.h"

namespace tudat
{

namespace propagators
{
//!Class for performing full numerical integration of a dynamical system in a single arc.
/*!
 *  Class for performing full numerical integration of a dynamical system in a single arc, i.e. the equations of motion
 *  have a single initial time, and are propagated once for the full prescribed time interval. This is in contrast to
 *  multi-arc dynamics, where the time interval si cut into pieces. In this class, the governing equations are set once,
 *  but can be re-integrated for different initial conditions using the same instance of the class.
 */
template< typename StateScalarType = double, typename TimeType = double >
class SingleArcDynamicsSimulator : public DynamicsSimulator< StateScalarType, TimeType >
{
public:
    SingleArcDynamicsSimulator(
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const bool areEquationsOfMotionToBeIntegrated = true,
            const PredefinedSingleArcStateDerivativeModels< StateScalarType, TimeType >& predefinedStateDerivativeModels =
                    PredefinedSingleArcStateDerivativeModels< StateScalarType, TimeType >( ),
            const bool isPartOfMultiArc = false );

    using DynamicsSimulator< StateScalarType, TimeType >::bodies_;

    //! Destructor
    ~SingleArcDynamicsSimulator( ) {}

    //! This function numerically (re-)integrates the equations of motion.
    /*!
     *  This function numerically (re-)integrates the equations of motion, using the settings set through the constructor
     *  and a new initial state vector provided here. The raw results are set in the equationsOfMotionNumericalSolution_
     *  \param initialStates Initial state vector that is to be used for numerical integration. Note that this state should
     *  be in the correct frame (i.e. corresponding to centralBodies in propagatorSettings_), but not in the propagator-
     *  specific form (i.e Encke, Gauss, etc. for translational dynamics)
     * \sa SingleStateTypeDerivative::convertToOutputSolution
     */
    void integrateEquationsOfMotion( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& initialStates );

    void integrate( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& initialStates );

    template< typename SimulationResults >
    void integrateEquationsOfMotion( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& processedInitialState,
                                     const std::shared_ptr< SimulationResults > propagationResults );

    //! Function to return the map of state history of numerically integrated bodies (base class interface).
    /*!
     * Function to return the map of state history of numerically integrated bodies (base class interface).
     * \return Vector is size 1, with entry: map of state history of numerically integrated bodies.
     */
    std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > getEquationsOfMotionNumericalSolutionBase( );

    //! Function to return the map of dependent variable history that was saved during numerical propagation (base class interface)
    /*!
     * Function to return the map of dependent variable history that was saved during numerical propagation (base class interface)
     * \return Vector is size 1, with entry: map of dependent variable history that was saved during numerical propagation.
     */
    std::vector< std::map< TimeType, Eigen::VectorXd > > getDependentVariableNumericalSolutionBase( );

    //! Function to return the map of cumulative computation time history that was saved during numerical propagation.
    /*!
     * Function to return the map of cumulative computation time history that was saved during numerical propagation (base class interface).
     * \return Vector is size 1, with entry: map of cumulative computation time history that was saved during numerical propagation.
     */
    std::vector< std::map< TimeType, double > > getCumulativeComputationTimeHistoryBase( );

    //! Function to get the settings for the numerical integrator.
    /*!
     * Function to get the settings for the numerical integrator.
     * \return The settings for the numerical integrator.
     */
    std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > getIntegratorSettings( )
    {
        return integratorSettings_;
    }

    //! Function to get the function that performs a single state derivative function evaluation.
    /*!
     * Function to get the function that performs a single state derivative function evaluation.
     * \return Function that performs a single state derivative function evaluation.
     */
    std::function< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >(
            const TimeType,
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& ) >
    getStateDerivativeFunction( )
    {
        return stateDerivativeFunction_;
    }

    //! Function to get the settings for the propagator.
    /*!
     * Function to get the settings for the propagator.
     * \return The settings for the propagator.
     */
    std::shared_ptr< SingleArcPropagatorSettings< StateScalarType, TimeType > > getPropagatorSettings( )
    {
        return propagatorSettings_;
    }

    //! Function to get the object that updates the environment.
    /*!
     * Function to get the object responsible for updating the environment based on the current state and time.
     * \return Object responsible for updating the environment based on the current state and time.
     */
    std::shared_ptr< EnvironmentUpdater< StateScalarType, TimeType > > getEnvironmentUpdater( )
    {
        return environmentUpdater_;
    }

    //! Function to get the object that updates and returns state derivative
    /*!
     * Function to get the object that updates current environment and returns state derivative from single function call
     * \return Object that updates current environment and returns state derivative from single function call
     */
    std::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > getDynamicsStateDerivative( )
    {
        return dynamicsStateDerivative_;
    }

    //! Function to retrieve the object defining when the propagation is to be terminated.
    /*!
     * Function to retrieve the object defining when the propagation is to be terminated.
     * \return Object defining when the propagation is to be terminated.
     */
    std::shared_ptr< PropagationTerminationCondition > getPropagationTerminationCondition( )
    {
        return propagationTerminationCondition_;
    }

    //! Function to retrieve the list of object that process the integrated numerical solution by updating the environment
    /*!
     * Function to retrieve the List of object (per dynamics type) that process the integrated numerical solution by
     * updating the environment
     * \return List of object (per dynamics type) that process the integrated numerical solution by updating the environment
     */
    std::map< IntegratedStateType, std::shared_ptr< SingleArcIntegratedStateProcessor< TimeType, StateScalarType > > >
    getIntegratedStateProcessors( )
    {
        return integratedStateProcessors_;
    }

    //! Function to retrieve initial time of propagation
    /*!
     * Function to retrieve initial time of propagation
     * \return Initial time of propagation
     */
    double getInitialPropagationTime( );

    //! Function to reset initial propagation time
    /*!
     * Function to reset initial propagation time
     * \param initialPropagationTime New initial propagation time
     */
    void resetInitialPropagationTime( const double initialPropagationTime );

    //! Function to retrieve the functions that compute the dependent variables at each time step
    /*!
     * Function to retrieve the functions that compute the dependent variables at each time step
     * \return Functions that compute the dependent variables at each time step
     */
    std::function< Eigen::VectorXd( ) > getDependentVariablesFunctions( )
    {
        return dependentVariablesFunctions_;
    }

    //! Function to reset the object that checks whether the simulation has finished from
    //! (newly defined) propagation settings.
    /*!
     *  Function to reset the object that checks whether the simulation has finished from
     *  (newly defined) propagation settings.
     */
    void resetPropagationTerminationConditions( );

    //! This function updates the environment with the numerical solution of the propagation.
    /*!
     *  This function updates the environment with the numerical solution of the propagation. It sets
     *  the propagated translational dynamics solution as the new input for the Ephemeris object of the body that was
     *  propagated.
     */
    void processNumericalEquationsOfMotionSolution( );

    void suppressDependentVariableDataPrinting( );

    void enableDependentVariableDataPrinting( );

    void createAndSetIntegratedStateProcessors( );

    std::shared_ptr< SimulationResults< StateScalarType, TimeType > > getPropagationResults( )
    {
        return propagationResults_;
    }

    std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType > > getSingleArcPropagationResults( )
    {
        return propagationResults_;
    }

    ///////////////////////////////////////////////////
    //////////////// DEPRECATED ///////////////////////
    ///////////////////////////////////////////////////

    //! Function to return the map of state history of numerically integrated bodies.
    /*!
     * Function to return the map of state history of numerically integrated bodies.
     * \return Map of state history of numerically integrated bodies.
     */
    const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& getEquationsOfMotionNumericalSolution( )
    {
        return propagationResults_->equationsOfMotionNumericalSolution_;
    }

    std::map< double, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > getEquationsOfMotionNumericalSolutionDouble( );

    //! Function to return the map of state history of numerically integrated bodies, in propagation coordinates.
    /*!
     * Function to return the map of state history of numerically integrated bodies, in propagation coordinates.
     * \return Map of state history of numerically integrated bodies, in propagation coordinates.
     */
    const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& getEquationsOfMotionNumericalSolutionRaw( )
    {
        return propagationResults_->equationsOfMotionNumericalSolutionRaw_;
    }

    std::map< double, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > getEquationsOfMotionNumericalSolutionRawDouble( );

    //! Function to return the map of dependent variable history that was saved during numerical propagation.
    /*!
     * Function to return the map of dependent variable history that was saved during numerical propagation.
     * \return Map of dependent variable history that was saved during numerical propagation.
     */
    const std::map< TimeType, Eigen::VectorXd >& getDependentVariableHistory( )
    {
        return propagationResults_->dependentVariableHistory_;
    }

    std::map< double, Eigen::VectorXd > getDependentVariableHistoryDouble( );

    //! Function to return the map of cumulative computation time history that was saved during numerical propagation.
    /*!
     * Function to return the map of cumulative computation time history that was saved during numerical propagation.
     * \return Map of cumulative computation time history that was saved during numerical propagation.
     */
    std::map< TimeType, double > getCumulativeComputationTimeHistory( )
    {
        return propagationResults_->cumulativeComputationTimeHistory_;
    }

    std::map< double, double > getCumulativeComputationTimeHistoryDouble( );

    //! Function to return the map of number of cumulative function evaluations that was saved during numerical propagation.
    /*!
     * Function to return the map of cumulative number of function evaluations that was saved during numerical propagation.
     * \return Map of cumulative number of function evaluations that was saved during numerical propagation.
     */
    std::map< TimeType, unsigned int > getCumulativeNumberOfFunctionEvaluations( )
    {
        return propagationResults_->cumulativeNumberOfFunctionEvaluations_;
    }

    std::map< double, unsigned int > getCumulativeNumberOfFunctionEvaluationsDouble( );

    //! Function to retrieve the event that triggered the termination of the last propagation
    /*!
     * Function to retrieve the event that triggered the termination of the last propagation
     * \return Event that triggered the termination of the last propagation
     */
    std::shared_ptr< PropagationTerminationDetails > getPropagationTerminationReason( )
    {
        return propagationResults_->propagationTerminationReason_;
    }

    void setPropagationTerminationReason( const std::shared_ptr< PropagationTerminationDetails > propagationTerminationReason )
    {
        propagationResults_->propagationTerminationReason_ = propagationTerminationReason;
    }

    //! Get whether the integration was completed successfully.
    /*!
     * Get whether the integration was completed successfully.
     * \return Whether the integration was completed successfully by reaching the termination condition.
     */
    virtual bool integrationCompletedSuccessfully( ) const;

    //! Function to retrieve the dependent variables IDs
    /*!
     * Function to retrieve the dependent variables IDs
     * \return Map listing starting entry of dependent variables in output vector, along with associated ID
     */
    std::map< std::pair< int, int >, std::string > getDependentVariableIds( );

    //! Function return whether the propagation is sequential or not (forward and backward leg).
    bool isPropagationSequential( ) const
    {
        return sequentialPropagation_;
    }

    ///////////////////////////////////////////////////
    //////////////// END DEPRECATED ///////////////////
    ///////////////////////////////////////////////////

    template< typename InputTimeType, typename InputStateScalarType >
    std::map< InputTimeType, Eigen::VectorXd > evaluateDependentVariablesAlongTrajectory(
            const std::map< InputTimeType, Eigen::Matrix< InputStateScalarType, Eigen::Dynamic, 1 > >& stateHistory );

protected:
    //! List of object (per dynamics type) that process the integrated numerical solution by updating the environment
    std::map< IntegratedStateType, std::shared_ptr< SingleArcIntegratedStateProcessor< TimeType, StateScalarType > > >
            integratedStateProcessors_;

    //! Object responsible for updating the environment based on the current state and time.
    /*!
     *  Object responsible for updating the environment based on the current state and time. Calling the updateEnvironment
     * function automatically updates all dependent variables that are needed to calulate the state derivative.
     */
    std::shared_ptr< EnvironmentUpdater< StateScalarType, TimeType > > environmentUpdater_;

    //! Interface object that updates current environment and returns state derivative from single function call.
    std::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > dynamicsStateDerivative_;

    //! Function that performs a single state derivative function evaluation.
    /*!
     *  Function that performs a single state derivative function evaluation, will typically be set to
     *  DynamicsStateDerivativeModel< TimeType, StateScalarType >::computeStateDerivative function.
     *  Calling this function will first update the environment (using environmentUpdater_) and then calculate the
     *  full system state derivative.
     */
    std::function< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >(
            const TimeType,
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& ) >
            stateDerivativeFunction_;

    //! Settings for numerical integrator.
    std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings_;

    //! Settings for propagator.
    std::shared_ptr< SingleArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings_;

    //! Object defining when the propagation is to be terminated.
    std::shared_ptr< PropagationTerminationCondition > propagationTerminationCondition_;

    //! Function returning dependent variables (during numerical propagation)
    std::function< Eigen::VectorXd( ) > dependentVariablesFunctions_;

    //    std::map< std::pair< int, int >, std::string > dependentVariableIds_;
    //
    //    std::map< std::pair< int, int >, std::string > processedStateIds_;
    //
    //    std::map< std::pair< int, int >, std::string > propagatedStateIds_;

    //! Object for retrieving ephemerides for transformation of reference frame (origins)
    std::shared_ptr< ephemerides::ReferenceFrameManager > frameManager_;

    std::shared_ptr< SingleArcPropagatorProcessingSettings > outputSettings_;

    std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType > > propagationResults_;

    //! Boolean denoting whether the propagation is performing sequentially, or both forward and backward (default = true).
    bool sequentialPropagation_;

private:
    //! Function that propagates the dynamics and (if requested) variational equations.
    /*
     *  Function that propagates the dynamics and (if requested) variational equations. Whether the variational
     *  equations are propagated is defined by the choice of SimulationResults template argument (if
     *  SingleArcSimulationResults< StateScalarType, TimeType >: dynamics only;
     *  if SingleArcVariationalSimulationResults< StateScalarType, TimeType >: dynamics and variational equations)
     *  NOTE: This function requires the performPropagationPreProcessingSteps
     *  and performPropagationPostProcessingSteps to be called before/after it. This is done automatically by the
     *  integrateEquationsOfMotion function.
     */
    template< typename SimulationResults >
    void propagateDynamics(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, SimulationResults::number_of_columns >& processedInitialState,
            const std::shared_ptr< SimulationResults > propagationResults );

    //! Function that propagates the dynamics and (if requested) variational equations until a given termination condition.
    /*
     *  Function that propagates the dynamics and (if requested) variational equations from the initial state until a given
     *  termination condition, and writes the results to the propagationResults object.
     */
    template< typename SimulationResults >
    void propagateUntilTermination(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, SimulationResults::number_of_columns >& processedInitialState,
            const std::shared_ptr< SimulationResults > propagationResults,
            const std::shared_ptr< PropagationTerminationCondition > terminationCondition );

    //! Function to perform steps necessary to reset all relevant models for the upcoming propagation
    /*
     *  Function to perform steps necessary to reset all relevant models for the upcoming propagation:
     *  - Whether to propagate dynamics and/or vatiational equations
     *  - Reset counter of function evaluations to zero
     *  - Reset termination conditions
     *  - Empty object holding the numerical simulation results of the previous run
     *  - Print messages to terminal, as requested by user settings
     */
    template< typename SimulationResults >
    void performPropagationPreProcessingSteps( const std::shared_ptr< SimulationResults > propagationResults );

    //! Function to perform steps necessary to finalize the propagation
    /*
     *  Function to perform steps necessary to finalize the propagation
     *  - Store number of function evaluations in the results object
     *  - Print messages to terminal, as requested by user settings
     *  - Update the environment (e.g. use numerical results to create tabulated ephemerides and similar for other dynamics)
     *    if requested by user
     */
    template< typename SimulationResults >
    void performPropagationPostProcessingSteps( const std::shared_ptr< SimulationResults > propagationResults );
};

#if TUDAT_BUILD_TESTS_WITH_EXPLICIT_INSTANTIATIONS
extern template class SingleArcDynamicsSimulator< double, double >;
#endif

}  // namespace propagators

}  // namespace tudat

#include "tudat/simulation/propagation_setup/singleArcDynamicsSimulatorImplementation.h"

#endif  // TUDAT_SINGLEARCDYNAMICSSIMULATOR_H
