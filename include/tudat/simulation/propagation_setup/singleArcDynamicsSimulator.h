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
            const std::shared_ptr< SingleArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings,
            const bool areEquationsOfMotionToBeIntegrated = true,
            const PredefinedSingleArcStateDerivativeModels< StateScalarType, TimeType >& predefinedStateDerivativeModels =
                    PredefinedSingleArcStateDerivativeModels< StateScalarType, TimeType >( ),
            const bool isPartOfMultiArc = false ):
        DynamicsSimulator< StateScalarType, TimeType >( bodies, propagatorSettings ), propagatorSettings_( propagatorSettings )
    {
        // Check consistency of input settings
        if( propagatorSettings == nullptr )
        {
            throw std::runtime_error( "Error in dynamics simulator, propagator settings not defined." );
        }
        else if( std::dynamic_pointer_cast< SingleArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ) == nullptr )
        {
            throw std::runtime_error( "Error in dynamics simulator, input must be single-arc." );
        }
        else
        {
            // Retrieve output and integrator settings TODO: no need to set as member variables; can just retrieve from propagatorSettings_
            outputSettings_ = propagatorSettings_->getOutputSettingsWithCheck( );
            integratorSettings_ = propagatorSettings_->getIntegratorSettings( );
        }
        if( integratorSettings_ == nullptr )
        {
            throw std::runtime_error( "Error in dynamics simulator, integrator settings not defined." );
        }
        checkPropagatedStatesFeasibility( propagatorSettings_, bodies_, isPartOfMultiArc );

        // Create objects that reset the environment (e.g. ephemerides) after propagation is required
        if( propagatorSettings_->getOutputSettings( )->getCreateStateProcessors( ) )
        {
            createAndSetIntegratedStateProcessors( );
        }

        // Create object that updates the environment during propagation
        try
        {
            environmentUpdater_ =
                    createEnvironmentUpdaterForDynamicalEquations< StateScalarType, TimeType >( propagatorSettings_, bodies_ );
        }
        catch( const std::runtime_error& error )
        {
            throw std::runtime_error( "Error when creating environment updater: " + std::string( error.what( ) ) );
        }

        // Create object that calculates the complete state derivatives
        if( predefinedStateDerivativeModels.stateDerivativeModels_.size( ) == 0 )
        {
            dynamicsStateDerivative_ = std::make_shared< DynamicsStateDerivativeModel< TimeType, StateScalarType > >(
                    createStateDerivativeModels< StateScalarType, TimeType >(
                            propagatorSettings_, bodies_, propagatorSettings_->getInitialTime( ) ),
                    std::bind( &EnvironmentUpdater< StateScalarType, TimeType >::updateEnvironment,
                               environmentUpdater_,
                               std::placeholders::_1,
                               std::placeholders::_2,
                               std::placeholders::_3 ) );
        }
        else
        {
            dynamicsStateDerivative_ = std::make_shared< DynamicsStateDerivativeModel< TimeType, StateScalarType > >(
                    predefinedStateDerivativeModels.stateDerivativeModels_,
                    std::bind( &EnvironmentUpdater< StateScalarType, TimeType >::updateEnvironment,
                               environmentUpdater_,
                               std::placeholders::_1,
                               std::placeholders::_2,
                               std::placeholders::_3 ) );
        }
        stateDerivativeFunction_ = std::bind( &DynamicsStateDerivativeModel< TimeType, StateScalarType >::computeStateDerivative,
                                              dynamicsStateDerivative_,
                                              std::placeholders::_1,
                                              std::placeholders::_2 );

        // Create object that determines if the propagation is to be terminated
        propagationTerminationCondition_ =
                createPropagationTerminationConditions( propagatorSettings_->getTerminationSettings( ),
                                                        bodies_,
                                                        integratorSettings_->initialTimeStep_,
                                                        dynamicsStateDerivative_->getStateDerivativeModels( ),
                                                        predefinedStateDerivativeModels.stateDerivativePartials_ );

        sequentialPropagation_ = true;
        if( propagationTerminationCondition_->getTerminationType( ) == non_sequential_stopping_condition )
        {
            sequentialPropagation_ = false;
            if( integratorSettings_->initialTimeStep_ < 0.0 )
            {
                throw std::runtime_error(
                        "Error when using non-sequential propagation, the initial integrator time step must be positive (first provided "
                        "for forward leg, "
                        "conversion to negative time step for backward leg is automatic)." );
            }
        }

        std::map< IntegratedStateType, std::vector< std::tuple< std::string, std::string, PropagatorType > > > integratedStateAndBodyList =
                getIntegratedTypeAndBodyList( propagatorSettings_ );

        std::map< std::pair< int, int >, std::string > dependentVariableIds_;
        std::map< std::pair< int, int >, std::shared_ptr< SingleDependentVariableSaveSettings > > orderedDependentVariableSettings_;

        // Create functions that compute the dependent variables
        if( propagatorSettings_->getDependentVariablesToSave( ).size( ) > 0 )
        {
            std::pair< std::function< Eigen::VectorXd( ) >, std::map< std::pair< int, int >, std::string > > dependentVariableData =
                    createDependentVariableListFunction< TimeType, StateScalarType >(
                            propagatorSettings_->getDependentVariablesToSave( ),
                            bodies_,
                            orderedDependentVariableSettings_,
                            dynamicsStateDerivative_->getStateDerivativeModels( ),
                            predefinedStateDerivativeModels.stateDerivativePartials_ );
            dependentVariablesFunctions_ = dependentVariableData.first;
            dependentVariableIds_ = dependentVariableData.second;
        }

        // Create object that will contain and process the propagation results
        std::shared_ptr< SingleArcDependentVariablesInterface< TimeType > > dependentVariableInterface =
                std::make_shared< SingleArcDependentVariablesInterface< TimeType > >(
                        std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, Eigen::VectorXd > >( ),
                        propagatorSettings_->getDependentVariablesToSave( ),
                        dependentVariableIds_,
                        orderedDependentVariableSettings_,
                        bodies );

        propagationResults_ = std::make_shared< SingleArcSimulationResults< StateScalarType, TimeType > >(
                integratedStateAndBodyList,
                propagatorSettings_->getOutputSettingsWithCheck( ),
                std::bind( &DynamicsStateDerivativeModel< TimeType, StateScalarType >::convertNumericalStateSolutionsToOutputSolutions,
                           dynamicsStateDerivative_,
                           std::placeholders::_1,
                           std::placeholders::_2 ),
                dependentVariableInterface,
                sequentialPropagation_ );

        // Integrate equations of motion if required.
        if( areEquationsOfMotionToBeIntegrated )
        {
            integrateEquationsOfMotion( propagatorSettings_->getInitialStates( ) );
        }
    }

    using DynamicsSimulator< StateScalarType, TimeType >::bodies_;

    //! Constructor of simulator.
    /*!
     *  Constructor of simulator, constructs integrator and object for calculating each time step of integration.
     *  \param bodies Map of bodies (with names) of all bodies in integration.
     *  \param integratorSettings Settings for numerical integrator.
     *  \param propagatorSettings Settings for propagator.
     *  \param areEquationsOfMotionToBeIntegrated Boolean to denote whether equations of motion should be integrated
     *  immediately at the end of the contructor or not (default true).
     *  \param clearNumericalSolutions Boolean to determine whether to clear the raw numerical solution member variables
     *  of this class, after propagation and resetting ephemerides (default false).
     *  \param setIntegratedResult Boolean to determine whether to automatically use the integrated results to set
     *  ephemerides (default false).
     *  \param printNumberOfFunctionEvaluations Boolean denoting whether the number of function evaluations should be printed
     *  at the end of propagation.
     *  \param initialClockTime Initial clock time from which to determine cumulative computation time.
     *  By default now( ), i.e. the moment at which this function is called.
     *  \param stateDerivativeModels List of state derivative models used in the simulation.
     */
    SingleArcDynamicsSimulator(
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const PredefinedSingleArcStateDerivativeModels< StateScalarType, TimeType >& predefinedStateDerivativeModels,
            const bool areEquationsOfMotionToBeIntegrated = true,
            const bool clearNumericalSolutions = false,
            const bool setIntegratedResult = false,
            const bool printNumberOfFunctionEvaluations = false,
            const std::chrono::steady_clock::time_point initialClockTime = std::chrono::steady_clock::now( ),
            const bool printDependentVariableData = false,
            const bool printStateData = false,
            const bool updateDependentVariableInterpolator = false ):
        SingleArcDynamicsSimulator( bodies,
                                    validateDeprecatedSingleArcSettings( integratorSettings,
                                                                         propagatorSettings,
                                                                         clearNumericalSolutions,
                                                                         setIntegratedResult,
                                                                         printNumberOfFunctionEvaluations,
                                                                         printDependentVariableData,
                                                                         printStateData,
                                                                         updateDependentVariableInterpolator ),
                                    areEquationsOfMotionToBeIntegrated,
                                    predefinedStateDerivativeModels )
    {}

    SingleArcDynamicsSimulator( const simulation_setup::SystemOfBodies& bodies,
                                const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
                                const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
                                const bool areEquationsOfMotionToBeIntegrated = true,
                                const bool clearNumericalSolutions = false,
                                const bool setIntegratedResult = false,
                                const bool printNumberOfFunctionEvaluations = false,
                                const bool printDependentVariableData = false,
                                const bool printStateData = false,
                                const bool updateDependentVariableInterpolator = false ):
        SingleArcDynamicsSimulator( bodies,
                                    integratorSettings,
                                    propagatorSettings,
                                    PredefinedSingleArcStateDerivativeModels< StateScalarType, TimeType >( ),
                                    areEquationsOfMotionToBeIntegrated,
                                    clearNumericalSolutions,
                                    setIntegratedResult,
                                    printNumberOfFunctionEvaluations,
                                    std::chrono::steady_clock::now( ),
                                    printDependentVariableData,
                                    printStateData,
                                    updateDependentVariableInterpolator )
    {}

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
    void integrateEquationsOfMotion( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& initialStates )
    {
        integrateEquationsOfMotion< SingleArcSimulationResults< StateScalarType, TimeType > >(
                dynamicsStateDerivative_->convertFromOutputSolution( initialStates, propagatorSettings_->getInitialTime( ) ),
                propagationResults_ );
    }

    void integrate( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& initialStates )
    {
        integrateEquationsOfMotion( initialStates );
    }

    template< typename SimulationResults >
    void integrateEquationsOfMotion( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& processedInitialState,
                                     const std::shared_ptr< SimulationResults > propagationResults )
    {
        performPropagationPreProcessingSteps( propagationResults );
        propagateDynamics< SimulationResults >(
                processedInitialState,
                propagationResults,
                PostProcessingFunctionProvider< StateScalarType, TimeType, SimulationResults::number_of_columns >::
                        getPostProcessingFunction( dynamicsStateDerivative_ ) );
        performPropagationPostProcessingSteps( propagationResults );
    }

    //! Function to return the map of state history of numerically integrated bodies (base class interface).
    /*!
     * Function to return the map of state history of numerically integrated bodies (base class interface).
     * \return Vector is size 1, with entry: map of state history of numerically integrated bodies.
     */
    std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > getEquationsOfMotionNumericalSolutionBase( )
    {
        return std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >(
                { getEquationsOfMotionNumericalSolution( ) } );
    }

    //! Function to return the map of dependent variable history that was saved during numerical propagation (base class interface)
    /*!
     * Function to return the map of dependent variable history that was saved during numerical propagation (base class interface)
     * \return Vector is size 1, with entry: map of dependent variable history that was saved during numerical propagation.
     */
    std::vector< std::map< TimeType, Eigen::VectorXd > > getDependentVariableNumericalSolutionBase( )
    {
        return std::vector< std::map< TimeType, Eigen::VectorXd > >( { getDependentVariableHistory( ) } );
    }

    //! Function to return the map of cumulative computation time history that was saved during numerical propagation.
    /*!
     * Function to return the map of cumulative computation time history that was saved during numerical propagation (base class interface).
     * \return Vector is size 1, with entry: map of cumulative computation time history that was saved during numerical propagation.
     */
    std::vector< std::map< TimeType, double > > getCumulativeComputationTimeHistoryBase( )
    {
        return std::vector< std::map< TimeType, double > >( { getCumulativeComputationTimeHistory( ) } );
    }

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
    double getInitialPropagationTime( )
    {
        return propagatorSettings_->getInitialTime( );
    }

    //! Function to reset initial propagation time
    /*!
     * Function to reset initial propagation time
     * \param initialPropagationTime New initial propagation time
     */
    void resetInitialPropagationTime( const double initialPropagationTime )

extern template class SingleArcDynamicsSimulator< double, double >;

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_SINGLEARCDYNAMICSSIMULATOR_H
