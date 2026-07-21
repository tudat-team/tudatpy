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
            std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > stateDerivativeModels =
                    createStateDerivativeModels< StateScalarType, TimeType >(
                            propagatorSettings_, bodies_, propagatorSettings_->getInitialTime( ) );
            stateDerivativeUpdater_ =
                    createStateDerivativeUpdaterForDynamicalEquations( propagatorSettings_, stateDerivativeModels, bodies_ );
            dynamicsStateDerivative_ = std::make_shared< DynamicsStateDerivativeModel< TimeType, StateScalarType > >(
                    stateDerivativeModels,
                    std::bind( &EnvironmentUpdater< StateScalarType, TimeType >::updateEnvironment,
                               environmentUpdater_,
                               std::placeholders::_1,
                               std::placeholders::_2,
                               std::placeholders::_3 ),
                    std::shared_ptr< VariationalEquations >( ),
                    stateDerivativeUpdater_ );
        }
        else
        {
            stateDerivativeUpdater_ = createStateDerivativeUpdaterForDynamicalEquations(
                    propagatorSettings_, predefinedStateDerivativeModels.stateDerivativeModels_, bodies_ );
            dynamicsStateDerivative_ = std::make_shared< DynamicsStateDerivativeModel< TimeType, StateScalarType > >(
                    predefinedStateDerivativeModels.stateDerivativeModels_,
                    std::bind( &EnvironmentUpdater< StateScalarType, TimeType >::updateEnvironment,
                               environmentUpdater_,
                               std::placeholders::_1,
                               std::placeholders::_2,
                               std::placeholders::_3 ),
                    std::shared_ptr< VariationalEquations >( ),
                    stateDerivativeUpdater_ );
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
                            predefinedStateDerivativeModels.stateDerivativePartials_,
                            predefinedStateDerivativeModels.parametersToEstimate_ );
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
    {
        propagatorSettings_->resetInitialTime( initialPropagationTime );
    }

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
    void resetPropagationTerminationConditions( )
    {
        propagationTerminationCondition_ = createPropagationTerminationConditions( propagatorSettings_->getTerminationSettings( ),
                                                                                   bodies_,
                                                                                   integratorSettings_->initialTimeStep_,
                                                                                   dynamicsStateDerivative_->getStateDerivativeModels( ) );
    }

    //! This function updates the environment with the numerical solution of the propagation.
    /*!
     *  This function updates the environment with the numerical solution of the propagation. It sets
     *  the propagated translational dynamics solution as the new input for the Ephemeris object of the body that was
     *  propagated.
     */
    void processNumericalEquationsOfMotionSolution( )
    {
        if( outputSettings_->getSetIntegratedResult( ) )
        {
            try
            {
                // Create and set interpolators for ephemerides
                resetIntegratedStates( propagationResults_->equationsOfMotionNumericalSolution_, integratedStateProcessors_ );
            }
            catch( const std::exception& caughtException )
            {
                std::cerr << "Error occured when post-processing single-arc integration results, and seting integrated states in "
                             "environment, caught error is: "
                          << std::endl
                          << std::endl;
                std::cerr << caughtException.what( ) << std::endl << std::endl;
                std::cerr << "The problem may be that there is an insufficient number of data points (epochs) at which propagation results "
                             "are produced. Integrated results are given at" +
                                std::to_string( propagationResults_->equationsOfMotionNumericalSolution_.size( ) ) + " epochs"
                          << std::endl;
            }

            // Clear numerical solution if so required.
            if( propagatorSettings_->getOutputSettings( )->getClearNumericalSolutions( ) )
            {
                propagationResults_->clearSolutionMaps( );
            }

            for( auto bodyIterator : bodies_.getMap( ) )
            {
                bodyIterator.second->updateConstantEphemerisDependentMemberQuantities( );
            }
        }
        else if( propagatorSettings_->getOutputSettings( )->getClearNumericalSolutions( ) )
        {
            propagationResults_->clearSolutionMaps( );
        }

        if( propagatorSettings_->getOutputSettings( )->getUpdateDependentVariableInterpolator( ) )
        {
            propagationResults_->updateDependentVariableInterface( );
        }
    }

    void suppressDependentVariableDataPrinting( )
    {
        outputSettings_->getPrintSettings( )->setPrintDependentVariableData( false );
    }

    void enableDependentVariableDataPrinting( )
    {
        outputSettings_->getPrintSettings( )->setPrintDependentVariableData( true );
    }

    void createAndSetIntegratedStateProcessors( )
    {
        frameManager_ = simulation_setup::createFrameManager( bodies_.getMap( ) );
        integratedStateProcessors_ =
                createIntegratedStateProcessors< TimeType, StateScalarType >( propagatorSettings_, bodies_, frameManager_ );
    }

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

    std::map< double, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > getEquationsOfMotionNumericalSolutionDouble( )
    {
        return propagationResults_->template getEquationsOfMotionNumericalSolutionTemplated< double >( );
    }

    //! Function to return the map of state history of numerically integrated bodies, in propagation coordinates.
    /*!
     * Function to return the map of state history of numerically integrated bodies, in propagation coordinates.
     * \return Map of state history of numerically integrated bodies, in propagation coordinates.
     */
    const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& getEquationsOfMotionNumericalSolutionRaw( )
    {
        return propagationResults_->equationsOfMotionNumericalSolutionRaw_;
    }

    std::map< double, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > getEquationsOfMotionNumericalSolutionRawDouble( )
    {
        return propagationResults_->template getEquationsOfMotionNumericalSolutionRawTemplated< double >( );
    }

    //! Function to return the map of dependent variable history that was saved during numerical propagation.
    /*!
     * Function to return the map of dependent variable history that was saved during numerical propagation.
     * \return Map of dependent variable history that was saved during numerical propagation.
     */
    const std::map< TimeType, Eigen::VectorXd >& getDependentVariableHistory( )
    {
        return propagationResults_->dependentVariableHistory_;
    }

    std::map< double, Eigen::VectorXd > getDependentVariableHistoryDouble( )
    {
        return propagationResults_->template getDependentVariableHistoryTemplated< double >( );
    }

    //! Function to return the map of cumulative computation time history that was saved during numerical propagation.
    /*!
     * Function to return the map of cumulative computation time history that was saved during numerical propagation.
     * \return Map of cumulative computation time history that was saved during numerical propagation.
     */
    std::map< TimeType, double > getCumulativeComputationTimeHistory( )
    {
        return propagationResults_->cumulativeComputationTimeHistory_;
    }

    std::map< double, double > getCumulativeComputationTimeHistoryDouble( )
    {
        return propagationResults_->getCumulativeComputationTimeHistory( );
    }

    //! Function to return the map of number of cumulative function evaluations that was saved during numerical propagation.
    /*!
     * Function to return the map of cumulative number of function evaluations that was saved during numerical propagation.
     * \return Map of cumulative number of function evaluations that was saved during numerical propagation.
     */
    std::map< TimeType, unsigned int > getCumulativeNumberOfFunctionEvaluations( )
    {
        return propagationResults_->cumulativeNumberOfFunctionEvaluations_;
    }

    std::map< double, unsigned int > getCumulativeNumberOfFunctionEvaluationsDouble( )
    {
        return propagationResults_->getCumulativeNumberOfFunctionEvaluations( );
    }

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
    virtual bool integrationCompletedSuccessfully( ) const
    {
        return ( propagationResults_->propagationTerminationReason_->getPropagationTerminationReason( ) == termination_condition_reached );
    }

    //! Function to retrieve the dependent variables IDs
    /*!
     * Function to retrieve the dependent variables IDs
     * \return Map listing starting entry of dependent variables in output vector, along with associated ID
     */
    std::map< std::pair< int, int >, std::string > getDependentVariableIds( )
    {
        return propagationResults_->getDependentVariableId( );
    }

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
            const std::map< InputTimeType, Eigen::Matrix< InputStateScalarType, Eigen::Dynamic, 1 > >& stateHistory )
    {
        std::map< InputTimeType, Eigen::VectorXd > dependentVariables_;

        // Integrate equations of motion numerically.
        simulation_setup::setAreBodiesInPropagation( bodies_, true );

        for( auto it : stateHistory )
        {
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > unprocessedState =
                    dynamicsStateDerivative_->convertFromOutputSolution( it.second.template cast< StateScalarType >( ), it.first );
            dynamicsStateDerivative_->computeStateDerivative( static_cast< double >( it.first ), unprocessedState );
            dependentVariables_[ it.first ] = dependentVariablesFunctions_( );
        }

        simulation_setup::setAreBodiesInPropagation( bodies_, false );
        return dependentVariables_;
    }

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

    std::shared_ptr< StateDerivativeUpdater< StateScalarType, TimeType > > stateDerivativeUpdater_;

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
            const std::shared_ptr< SimulationResults > propagationResults,
            const std::function< void( Eigen::Matrix< StateScalarType, Eigen::Dynamic, SimulationResults::number_of_columns >& ) >
                    statePostProcessingFunction )
    {
        // Integrate equations of motion numerically.
        simulation_setup::setAreBodiesInPropagation( bodies_, true );
        dynamicsStateDerivative_->updateStateDerivativeModelSettings(
                processedInitialState.block( 0, processedInitialState.cols( ) - 1, processedInitialState.rows( ), 1 ) );

        if( sequentialPropagation_ )
        {
            integrateEquations< SimulationResults,
                                Eigen::Matrix< StateScalarType, Eigen::Dynamic, SimulationResults::number_of_columns >,
                                TimeType >( stateDerivativeFunction_,
                                            processedInitialState,
                                            propagatorSettings_->getInitialTime( ),
                                            integratorSettings_,
                                            propagationTerminationCondition_,
                                            propagationResults,
                                            dependentVariablesFunctions_,
                                            statePostProcessingFunction,
                                            propagatorSettings_->getOutputSettings( ),
                                            dynamicsStateDerivative_ );
        }
        else
        {
            std::shared_ptr< NonSequentialPropagationTerminationCondition > nonSequentialTerminations =
                    std::dynamic_pointer_cast< NonSequentialPropagationTerminationCondition >( propagationTerminationCondition_ );
            integrateEquations< SimulationResults,
                                Eigen::Matrix< StateScalarType, Eigen::Dynamic, SimulationResults::number_of_columns >,
                                TimeType >( stateDerivativeFunction_,
                                            processedInitialState,
                                            propagatorSettings_->getInitialTime( ),
                                            integratorSettings_,
                                            nonSequentialTerminations->getForwardPropagationTerminationCondition( ),
                                            propagationResults,
                                            dependentVariablesFunctions_,
                                            statePostProcessingFunction,
                                            propagatorSettings_->getOutputSettings( ),
                                            dynamicsStateDerivative_ );

            integratorSettings_->initialTimeStep_ *= -1.0;
            integrateEquations< SimulationResults,
                                Eigen::Matrix< StateScalarType, Eigen::Dynamic, SimulationResults::number_of_columns >,
                                TimeType >( stateDerivativeFunction_,
                                            processedInitialState,
                                            propagatorSettings_->getInitialTime( ),
                                            integratorSettings_,
                                            nonSequentialTerminations->getBackwardPropagationTerminationCondition( ),
                                            propagationResults,
                                            dependentVariablesFunctions_,
                                            statePostProcessingFunction,
                                            propagatorSettings_->getOutputSettings( ),
                                            dynamicsStateDerivative_ );
            integratorSettings_->initialTimeStep_ *= -1.0;
        }

        simulation_setup::setAreBodiesInPropagation( bodies_, false );
    }

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
    void performPropagationPreProcessingSteps( const std::shared_ptr< SimulationResults > propagationResults )
    {
        // Reset functions
        dynamicsStateDerivative_->setPropagationSettings( std::vector< IntegratedStateType >( ), true, SimulationResults::is_variational );
        dynamicsStateDerivative_->resetFunctionEvaluationCounter( );
        dynamicsStateDerivative_->resetCumulativeFunctionEvaluationCounter( );
        resetPropagationTerminationConditions( );

        // Empty solution maps
        propagationResults->reset( );

        PropagationPrintingInterface< SimulationResults, StateScalarType, TimeType >::printSingleArcPrePropagationMessages(
                outputSettings_->getPrintSettings( ), outputSettings_->getPropagationStartHeader( ), propagationResults );
    }

    //! Function to perform steps necessary to finalize the propagation
    /*
     *  Function to perform steps necessary to finalize the propagation
     *  - Store number of function evaluations in the results object
     *  - Print messages to terminal, as requested by user settings
     *  - Update the environment (e.g. use numerical results to create tabulated ephemerides and similar for other dynamics)
     *    if requested by user
     */
    template< typename SimulationResults >
    void performPropagationPostProcessingSteps( const std::shared_ptr< SimulationResults > propagationResults )
    {
        // Retrieve number of cumulative function evaluations
        propagationResults->finalizePropagation( dynamicsStateDerivative_->getCumulativeNumberOfFunctionEvaluations( ) );
        PropagationPrintingInterface< SimulationResults, StateScalarType, TimeType >::printSingleArcPostPropagationMessages(
                outputSettings_->getPrintSettings( ), outputSettings_->getPropagationEndHeader( ), propagationResults );
        processNumericalEquationsOfMotionSolution( );
    }
};

#if TUDAT_BUILD_EXPLICIT_INSTANTIATIONS
extern template class SingleArcDynamicsSimulator< double, double >;
#endif

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_SINGLEARCDYNAMICSSIMULATOR_H
