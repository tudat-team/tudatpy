/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_HYBRIDARCVARIATIONALEQUATIONSSOLVER_H
#define TUDAT_HYBRIDARCVARIATIONALEQUATIONSSOLVER_H

#include "tudat/simulation/estimation_setup/singleArcVariationalEquationsSolver.h"
#include "tudat/simulation/estimation_setup/multiArcVariationalEquationsSolver.h"

namespace tudat
{

namespace propagators
{
//! Class to manage and execute the numerical integration of variational equations of a dynamical system in a combination
//! of single and multiple arcs
/*!
 *  Class to manage and execute the numerical integration of variational equations of a dynamical system, in addition
 *  to the dynamics itself, in a combination  of single and multiple arcs. In this class, the governing equations are set once,
 *  but can be re-integrated for different initial conditions using the same instance of the class.
 */
template< typename StateScalarType = double, typename TimeType = double >
class HybridArcVariationalEquationsSolver : public VariationalEquationsSolver< StateScalarType, TimeType >
{
public:
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > MatrixType;
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > VectorType;
    typedef HybridArcSimulationResults< SingleArcVariationalSimulationResults, StateScalarType, TimeType > HybridArcResults;

    using VariationalEquationsSolver< StateScalarType, TimeType >::parametersToEstimate_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::bodies_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::stateTransitionMatrixSize_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::parameterVectorSize_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::stateTransitionInterface_;

    HybridArcVariationalEquationsSolver(
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< HybridArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
            const bool integrateEquationsOnCreation = false ):
        VariationalEquationsSolver< StateScalarType, TimeType >(
                bodies,
                parametersToEstimate,
                propagatorSettings != nullptr ? propagatorSettings->getOutputSettingsWithCheck( )->getClearNumericalSolutions( ) : false )
    {
        initializeHybridArcVariationalEquationsSolver( bodies, propagatorSettings, integrateEquationsOnCreation );
    }

    //! Constructor
    /*!
     *  Constructor, sets up object for automatic evaluation and numerical integration of variational equations and equations of motion.
     *  \param bodies Map of bodies (with names) of all bodies in integration.
     *  \param integratorSettings Settings for numerical integrator.
     *  \param propagatorSettings Settings for propagator.
     *  \param parametersToEstimate Object containing all parameters that are to be estimated and their current settings and values.
     *  \param arcStartTimes Start times for separate arcs
     *  \param integrateDynamicalAndVariationalEquationsConcurrently Boolean defining whether variational and dynamical
     *  equations are to be propagated concurrently (if true) or sequentially (of false)
     *  \param clearNumericalSolution Boolean to determine whether to clear the raw numerical solution member variables
     *  (default true) after propagation and resetting of state transition interface.
     *  \param integrateEquationsOnCreation Boolean to denote whether equations should be integrated immediately at the
     *  end of this contructor.
     */
    HybridArcVariationalEquationsSolver(
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
            const std::vector< double > arcStartTimes,
            const bool integrateDynamicalAndVariationalEquationsConcurrently = true,
            const bool clearNumericalSolution = false,
            const bool integrateEquationsOnCreation = false,
            const bool setDependentVariablesInterface = false ):
        HybridArcVariationalEquationsSolver< StateScalarType, TimeType >(
                bodies,
                validateDeprecatedHybridArcSettings< StateScalarType, TimeType >( integratorSettings,
                                                                                  propagatorSettings,
                                                                                  arcStartTimes,
                                                                                  clearNumericalSolution,
                                                                                  true,
                                                                                  setDependentVariablesInterface ),
                parametersToEstimate,
                integrateEquationsOnCreation )
    {}

    HybridArcVariationalEquationsSolver(
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > singleArcIntegratorSettings,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > multiArcIntegratorSettings,
            const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
            const std::vector< double > arcStartTimes,
            const bool integrateDynamicalAndVariationalEquationsConcurrently = true,
            const bool clearNumericalSolution = false,
            const bool integrateEquationsOnCreation = false,
            const bool setDependentVariablesInterface = false ):
        HybridArcVariationalEquationsSolver< StateScalarType, TimeType >(
                bodies,
                validateDeprecatedHybridArcSettings< StateScalarType, TimeType >( singleArcIntegratorSettings,
                                                                                  multiArcIntegratorSettings,
                                                                                  propagatorSettings,
                                                                                  arcStartTimes,
                                                                                  clearNumericalSolution,
                                                                                  true,
                                                                                  setDependentVariablesInterface ),
                parametersToEstimate,
                integrateEquationsOnCreation )
    {}

    void initializeHybridArcVariationalEquationsSolver( const simulation_setup::SystemOfBodies& bodies,
                                                        const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
                                                        const bool integrateEquationsOnCreation )
    {
        // Cast propagator settings to correct type and check validity
        originalPopagatorSettings_ =
                std::dynamic_pointer_cast< HybridArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings );
        if( originalPopagatorSettings_ == nullptr )
        {
            throw std::runtime_error(
                    "Error when making HybridArcVariationalEquationsSolver, input propagation settings are not hybrid arc" );
        }

        // Retrive arc properties
        singleArcInitialTime_ = originalPopagatorSettings_->getSingleArcPropagatorSettings( )->getInitialTime( );
        int numberOfArcs = originalPopagatorSettings_->getMultiArcPropagatorSettings( )->getNmberOfArcs( );
        arcStartTimes_ = estimatable_parameters::getMultiArcStateEstimationArcStartTimes( parametersToEstimate_, false );

        // Get input size of single-arc and input multi-arc
        singleArcDynamicsSize_ = originalPopagatorSettings_->getSingleArcPropagatorSettings( )->getConventionalStateSize( );
        originalMultiArcDynamicsSize_ = originalPopagatorSettings_->getMultiArcPropagatorSettings( )->getConventionalStateSize( );
        for( unsigned int i = 0; i < arcStartTimes_.size( ); i++ )
        {
            originalMultiArcDynamicsSingleArcSize_.push_back( originalPopagatorSettings_->getMultiArcPropagatorSettings( )
                                                                      ->getSingleArcSettings( )
                                                                      .at( i )
                                                                      ->getConventionalStateSize( ) );
        }

        // Create propagator settings with the single arc settings included (at the beginning) in each arc
        std::shared_ptr< MultiArcPropagatorSettings< StateScalarType, TimeType > > extendedMultiArcSettings =
                getExtendedMultiPropagatorSettings( originalPopagatorSettings_->getSingleArcPropagatorSettings( ),
                                                    originalPopagatorSettings_->getMultiArcPropagatorSettings( ),
                                                    numberOfArcs );

        multiArcDynamicsSize_ = extendedMultiArcSettings->getConventionalStateSize( );
        for( unsigned int i = 0; i < arcStartTimes_.size( ); i++ )
        {
            multiArcDynamicsSingleArcSize_.push_back(
                    extendedMultiArcSettings->getSingleArcSettings( ).at( i )->getConventionalStateSize( ) );
        }

        propagatorSettings_ = std::make_shared< HybridArcPropagatorSettings< StateScalarType, TimeType > >(
                originalPopagatorSettings_->getSingleArcPropagatorSettings( )->clone( ), extendedMultiArcSettings );

        // Update estimated parameter vector to extended multi-arc settings
        setExtendedMultiArcParameters( arcStartTimes_ );

        // Create multi-arc solver with original parameter set
        extendedMultiArcSettings->getOutputSettings( )->setClearNumericalSolutions( false );
        extendedMultiArcSettings->getOutputSettings( )->setIntegratedResult( false );

        originalMultiArcSolver_ = std::make_shared< MultiArcVariationalEquationsSolver< StateScalarType, TimeType > >(
                bodies, originalPopagatorSettings_->getMultiArcPropagatorSettings( ), originalMultiArcParametersToEstimate_, false );

        // Create variational equations solvers for single- and multi-arc
        singleArcSolver_ = std::make_shared< SingleArcVariationalEquationsSolver< StateScalarType, TimeType > >(
                bodies, originalPopagatorSettings_->getSingleArcPropagatorSettings( ), singleArcParametersToEstimate_, true, false );

        multiArcSolver_ = std::make_shared< MultiArcVariationalEquationsSolver< StateScalarType, TimeType > >(
                bodies, extendedMultiArcSettings, multiArcParametersToEstimate_, false );

        for( unsigned int i = 0; i < multiArcSolver_->getDynamicsStateDerivatives( ).size( ); i++ )
        {
            multiArcSolver_->getDynamicsStateDerivatives( ).at( i )->getVariationalEquationsCalculator( )->suppressParameterCoupling(
                    propagatorSettings_->getSingleArcPropagatorSettings( )->getPropagatedStateSize( ) );
        }

        // Create function to retrieve single-arc initial states for extended multi-arc
        std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType, TimeType > > singleArcPropagationSettings =
                std::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType, TimeType > >(
                        propagatorSettings_->getSingleArcPropagatorSettings( ) );
        if( singleArcPropagationSettings == nullptr )
        {
            throw std::runtime_error( "Error when making HybridArcVariationalEquationsSolver, input single arc is not translational" );
        }

        variationalPropagationResults_ =
                std::make_shared< HybridArcResults >( singleArcSolver_->getSingleArcVariationalPropagationResults( ),
                                                      originalMultiArcSolver_->getMultiArcVariationalPropagationResults( ) );
        initialStatesFromSingleArcPropagation_ = std::bind( &getInitialStatesOfBodiesFromFrameManager< TimeType, StateScalarType >,
                                                            singleArcPropagationSettings->bodiesToIntegrate_,
                                                            singleArcPropagationSettings->centralBodies_,
                                                            bodies,
                                                            std::placeholders::_1,
                                                            createFrameManager( bodies.getMap( ) ) );

        // Propagate dynamical equations if requested
        if( integrateEquationsOnCreation )
        {
            integrateVariationalAndDynamicalEquations( propagatorSettings_->getInitialStates( ), 1 );
        }
    }

    //! Destructor
    ~HybridArcVariationalEquationsSolver( ) {}

    //! Function to integrate variational equations and equations of motion.
    /*!
     *  Function to integrate variational equations and equations of motion. At the end of this function,
     *  the stateTransitionInterface_ is reset with the new state transition and sensitivity matrices. If dynamical
     *  solution is to be processed, the environment is also updated to the new solution.
     *  \param initialStateEstimate Initial statez of the equations of motion that is to be used (in same order as in
     *  parametersToEstimate_). The initial states of single and multi-arcs propagations are concatenated into a single
     *  vector.
     *  \param integrateEquationsConcurrently Variable determining whether the equations of motion are to be
     *  propagated concurrently with variational equations of motion (if true), or before variational equations (if false).
     */
    void integrateVariationalAndDynamicalEquations( const VectorType& initialStateEstimate, const bool integrateEquationsConcurrently )
    {
        // TODO: do process depdendent variables in original multi-arc solver, do not process dependent variables in
        // extended solver. Also add dependent variables to original multi-arc solver.

        // Reset initial time and propagate multi-arc equations
        singleArcSolver_->integrateVariationalAndDynamicalEquations( initialStateEstimate.block( 0, 0, singleArcDynamicsSize_, 1 ),
                                                                     integrateEquationsConcurrently );

        // Extract single arc state to update multi-arc initial states
        resetMultiArcInitialStates( initialStateEstimate.block( singleArcDynamicsSize_, 0, multiArcDynamicsSize_, 1 ) );

        // Reset initial time and propagate single-arc equations
        multiArcSolver_->integrateVariationalAndDynamicalEquations(
                propagatorSettings_->getMultiArcPropagatorSettings( )->getInitialStates( ), integrateEquationsConcurrently );

        copyExtendedMultiArcInitialStatesToOriginalSettins( );

        // Extract multi-arc solution of dynamics, and remove the single arc bodies from the map.
        std::vector< std::map< TimeType, VectorType > > numericalMultiArcSolution =
                multiArcSolver_->getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );
        std::vector< std::map< TimeType, Eigen::VectorXd > > dependentVariableHistory =
                multiArcSolver_->getDynamicsSimulator( )->getDependentVariableHistory( );

        removeSingleArcBodiesFromMultiArcSolultion( numericalMultiArcSolution );

        originalMultiArcSolver_->getDynamicsSimulator( )->getMultiArcPropagationResults( )->restartPropagation( );

        // Reset original multi-arc bodies' dynamics
        originalMultiArcSolver_->getDynamicsSimulator( )->getMultiArcPropagationResults( )->manuallySetPropagationResults(
                numericalMultiArcSolution );
        originalMultiArcSolver_->getDynamicsSimulator( )->getMultiArcPropagationResults( )->manuallySetSecondaryData(
                multiArcSolver_->getDynamicsSimulator( )->getMultiArcPropagationResults( ) );

        originalMultiArcSolver_->getDynamicsSimulator( )->processNumericalEquationsOfMotionSolution( );

        // Create state transition matrix if not yet created.
        if( stateTransitionInterface_ == nullptr )
        {
            if( std::dynamic_pointer_cast< SingleArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                        singleArcSolver_->getStateTransitionMatrixInterface( ) ) == nullptr )
            {
                throw std::runtime_error( "Error when making hybrid state transition/sensitivity interface, single-arc input is nullptr" );
            }

            if( std::dynamic_pointer_cast< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface< StateScalarType > >(
                        multiArcSolver_->getStateTransitionMatrixInterface( ) ) == nullptr )
            {
                throw std::runtime_error( "Error when making hybrid state transition/sensitivity interface, multi-arc input is nullptr" );
            }

            stateTransitionInterface_ =
                    std::make_shared< HybridArcCombinedStateTransitionAndSensitivityMatrixInterface< StateScalarType > >(
                            std::dynamic_pointer_cast< SingleArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                                    singleArcSolver_->getStateTransitionMatrixInterface( ) ),
                            std::dynamic_pointer_cast< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface< StateScalarType > >(
                                    multiArcSolver_->getStateTransitionMatrixInterface( ) ) );
        }
    }

    //! Function to integrate equations of motion only.
    /*!
     *  Function to integrate equations of motion only.  If dynamical
     *  solution is to be processed, the environment is also updated to the new solution.
     *  \param initialStateEstimate Initial state of the equations of motion that is to be used (in same order as in
     *  parametersToEstimate_). The initial states of single and multi-arcs propagations are concatenated into a single
     *  vector.
     */
    void integrateDynamicalEquationsOfMotionOnly( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialStateEstimate )
    {
        // Reset initial time and propagate multi-arc equations
        singleArcSolver_->integrateDynamicalEquationsOfMotionOnly( initialStateEstimate.block( 0, 0, singleArcDynamicsSize_, 1 ) );

        // Extract single arc state to update multi-arc initial states
        resetMultiArcInitialStates( initialStateEstimate.block( singleArcDynamicsSize_, 0, multiArcDynamicsSize_, 1 ) );

        // Reset initial time and propagate single-arc equations
        multiArcSolver_->integrateDynamicalEquationsOfMotionOnly(
                propagatorSettings_->getMultiArcPropagatorSettings( )->getInitialStates( ) );

        copyExtendedMultiArcInitialStatesToOriginalSettins( );

        // Extract multi-arc solution of dynamics, and remove the single arc bodies from the map.
        std::vector< std::map< TimeType, VectorType > > numericalMultiArcSolution =
                multiArcSolver_->getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );
        std::vector< std::map< TimeType, Eigen::VectorXd > > dependentVariableHistory =
                multiArcSolver_->getDynamicsSimulator( )->getDependentVariableHistory( );

        removeSingleArcBodiesFromMultiArcSolultion( numericalMultiArcSolution );

        // Reset original multi-arc bodies' dynamics
        originalMultiArcSolver_->getDynamicsSimulator( )->getMultiArcPropagationResults( )->manuallySetPropagationResults(
                numericalMultiArcSolution );
        originalMultiArcSolver_->getDynamicsSimulator( )->getMultiArcPropagationResults( )->manuallySetSecondaryData(
                multiArcSolver_->getDynamicsSimulator( )->getMultiArcPropagationResults( ) );
        originalMultiArcSolver_->getDynamicsSimulator( )->processNumericalEquationsOfMotionSolution( );
    }

    //! Function to reset parameter estimate and re-integrate equations of motion and, if desired, variational equations.
    /*!
     *  Function to reset parameter estimate and re-integrate equations of motion and, if desired, variational equations
     *  using the new physical parameters/body initial states.
     *  \param newParameterEstimate New estimate of parameters that are to be estimated, in same order as defined
     *  in parametersToEstimate_ member.
     *  \param areVariationalEquationsToBeIntegrated Boolean defining whether the variational equations are to be
     *  reintegrated with the new parameter values.
     */
    void resetParameterEstimate( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > newParameterEstimate,
                                 const bool areVariationalEquationsToBeIntegrated = true )
    {
        // Reset values of parameters.
        parametersToEstimate_->template resetParameterValues< StateScalarType >( newParameterEstimate );
        simulation_setup::setInitialStateVectorFromParameterSet< StateScalarType >( parametersToEstimate_, originalPopagatorSettings_ );

        propagatorSettings_->getSingleArcPropagatorSettings( )->resetInitialStates(
                newParameterEstimate.segment( 0, singleArcDynamicsSize_ ) );
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > totalMultiArcInitialState =
                propagatorSettings_->getMultiArcPropagatorSettings( )->getInitialStates( );

        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > newInitialStatesOriginalPropagatorSettings =
                originalPopagatorSettings_->getInitialStates( );

        unsigned int counterFullArcWiseIndex = 0;
        unsigned int counterOriginalArcWiseIndex = 0;
        for( unsigned int i = 0; i < propagatorSettings_->getMultiArcPropagatorSettings( )->getSingleArcSettings( ).size( ); i++ )
        {
            totalMultiArcInitialState.segment( counterFullArcWiseIndex, singleArcDynamicsSize_ ) =
                    TUDAT_NAN * Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Ones( singleArcDynamicsSize_ );
            totalMultiArcInitialState.segment(
                    counterFullArcWiseIndex /*i * multiArcDynamicsSingleArcSize_.at( i )*/ + singleArcDynamicsSize_,
                    originalMultiArcDynamicsSingleArcSize_.at( i ) ) =
                    newInitialStatesOriginalPropagatorSettings.segment(
                            singleArcDynamicsSize_ + counterOriginalArcWiseIndex /*i * originalMultiArcDynamicsSingleArcSize_.at( i )*/,
                            originalMultiArcDynamicsSingleArcSize_.at( i ) );

            counterFullArcWiseIndex += multiArcDynamicsSingleArcSize_.at( i );
            counterOriginalArcWiseIndex += originalMultiArcDynamicsSingleArcSize_.at( i );
        }
        propagatorSettings_->getMultiArcPropagatorSettings( )->resetInitialStates( totalMultiArcInitialState );
        propagatorSettings_->updateInitialState( );

        // Reset parameters for arc-wise parameters in both originalMultiArcSolver_ and multiArcSolver_
        for( unsigned int i = 0; i < arcStartTimes_.size( ); i++ )
        {
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > newParametersValues =
                    originalPopagatorSettings_->getMultiArcPropagatorSettings( )->getSingleArcSettings( ).at( i )->getInitialStates( );
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > newArcWiseParametersValues =
                    originalMultiArcSolver_->getArcWiseParametersToEstimate( )
                            .at( i )
                            ->template getFullParameterValues< StateScalarType >( );
            newArcWiseParametersValues.segment( 0, newParametersValues.size( ) ) = newParametersValues;
            originalMultiArcSolver_->getArcWiseParametersToEstimate( ).at( i )->resetParameterValues( newArcWiseParametersValues );

            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > newFullParametersValues =
                    propagatorSettings_->getMultiArcPropagatorSettings( )->getSingleArcSettings( ).at( i )->getInitialStates( );
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > newFullArcWiseParametersValues =
                    multiArcSolver_->getArcWiseParametersToEstimate( ).at( i )->template getFullParameterValues< StateScalarType >( );
            newFullArcWiseParametersValues.segment( 0, newFullParametersValues.size( ) ) = newFullParametersValues;
            multiArcSolver_->getArcWiseParametersToEstimate( ).at( i )->resetParameterValues( newFullArcWiseParametersValues );
        }

        // Check if re-integration of variational equations is requested
        if( areVariationalEquationsToBeIntegrated )
        {
            // Integrate variational and state equations.
            this->integrateVariationalAndDynamicalEquations( propagatorSettings_->getInitialStates( ), 1 );
        }
        else
        {
            this->integrateDynamicalEquationsOfMotionOnly( propagatorSettings_->getInitialStates( ) );
        }
    }

    //! Function to retrieve propagator settings used for equations of motion
    /*!
     * Function to retrieve propagator settings used for equations of motion
     * \return Propagator settings used for equations of motion
     */
    std::shared_ptr< HybridArcPropagatorSettings< StateScalarType, TimeType > > getPropagatorSettings( )
    {
        return propagatorSettings_;
    }

    //! Function to retrieve the dynamics simulator object (as base-class pointer)
    /*!
     * Function to retrieve the dynamics simulator object (as base-class pointer). This function is not yet implemented
     * in hybric-arc model, as no single DynamicsSimulator model is used. Calling this function throws an error
     * \return Dynamics simulator object (as base-class pointer)
     */
    virtual std::shared_ptr< DynamicsSimulator< StateScalarType, TimeType > > getDynamicsSimulatorBase( )
    {
        throw std::runtime_error( "Error, getDynamicsSimulatorBase not implemented in hyrbid arc propagator" );
    }

    std::shared_ptr< MultiArcVariationalEquationsSolver< StateScalarType, TimeType > > getMultiArcSolver( )
    {
        return multiArcSolver_;
    }

    //! Object to solve single-arc variational equations.
    std::shared_ptr< SingleArcVariationalEquationsSolver< StateScalarType, TimeType > > getSingleArcSolver( )
    {
        return singleArcSolver_;
    }

    std::shared_ptr< HybridArcResults > getHybridArcVariationalPropagationResults( )
    {
        return variationalPropagationResults_;
    }

    std::shared_ptr< SimulationResults< StateScalarType, TimeType > > getVariationalPropagationResults( )
    {
        return getHybridArcVariationalPropagationResults( );
    }

protected:
    //! Function to set and process the arc start times of the multi-arc propagation
    /*!
     * Function to set and process the arc start times of the multi-arc propagation
     * \param arcStartTimes Arc start times of the multi-arc propagation
     */
    void setExtendedMultiArcParameters( const std::vector< double >& arcStartTimes )
    {
        // Retrieve and set original single and multi-arc parameter set
        singleArcParametersToEstimate_ = createEstimatableParameterSetArcSubSet( parametersToEstimate_, true );
        originalMultiArcParametersToEstimate_ = createEstimatableParameterSetArcSubSet( parametersToEstimate_, false );

        std::vector<
                std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > >
                singleArcParameters = parametersToEstimate_->getEstimatedSingleArcInitialStateParameters( );
        std::vector<
                std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > >
                originalMultiArcParameters = parametersToEstimate_->getEstimatedMultiArcInitialStateParameters( );

        // Get multi-arc parameters associated with estimated single-arc parameters
        std::vector<
                std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > >
                extendedMultiArcParameters;
        for( unsigned int i = 0; i < singleArcParameters.size( ); i++ )
        {
            extendedMultiArcParameters.push_back(
                    simulation_setup::getAssociatedMultiArcParameter( singleArcParameters.at( i ), arcStartTimes ) );
        }

        // Add original multi-arc parameters
        for( unsigned int i = 0; i < originalMultiArcParameters.size( ); i++ )
        {
            extendedMultiArcParameters.push_back( originalMultiArcParameters.at( i ) );
        }

        // Create multi-arc parameter set with single-arc parameters extended into multi-arc
        multiArcParametersToEstimate_ = std::make_shared< estimatable_parameters::EstimatableParameterSet< StateScalarType > >(
                parametersToEstimate_->getEstimatedDoubleParameters( ),
                parametersToEstimate_->getEstimatedVectorParameters( ),
                extendedMultiArcParameters );
        //        std::cout << "TEST: " << "\n\n";
        estimatable_parameters::printEstimatableParameterEntries( multiArcParametersToEstimate_ );
    }

    //! Function to reset the initial multi-arc states
    /*!
     * Function to reset the initial multi-arc states
     * \param manualMultiArcStates New multi-arc states
     */
    void resetMultiArcInitialStates( const VectorType& manualMultiArcStates )
    {
        //        std::cout << "manualMultiArcStates: " << manualMultiArcStates.transpose( ) << "\n\n";
        //        std::cout << "size manualMultiArcStates: " << manualMultiArcStates.size( ) << "\n\n";

        // Retrieve full multi-arc initial states, with single-arc bodies not (correctly) set
        std::vector< VectorType > arcInitialStates = propagatorSettings_->getMultiArcPropagatorSettings( )->getInitialStateList( );

        // Retrieve single-arc states from ephemerides
        int currentArcSize = 0;
        int indexMultiArcState = 0;
        for( unsigned int i = 0; i < arcInitialStates.size( ); i++ )
        {
            currentArcSize = arcInitialStates[ i ].rows( );
            //            std::cout << "currentArcSize: " << currentArcSize << "\n\n";
            //            std::cout << "singleArcDynamicsSize_" << singleArcDynamicsSize_ << "\n\n";
            arcInitialStates[ i ].segment( 0, singleArcDynamicsSize_ ) = initialStatesFromSingleArcPropagation_( arcStartTimes_.at( i ) );
            //            std::cout << "new single arc initial state:" << initialStatesFromSingleArcPropagation_( arcStartTimes_.at( i )
            //            ).transpose( ) << "\n\n";

            //            std::cout << "currentArcSize - singleArcDynamicsSize_: " << currentArcSize - singleArcDynamicsSize_ << "\n\n";
            //            std::cout << "i * currentArcSize + singleArcDynamicsSize_: " << i * currentArcSize + singleArcDynamicsSize_ <<
            //            "\n\n";
            arcInitialStates[ i ].segment( singleArcDynamicsSize_, currentArcSize - singleArcDynamicsSize_ ) = manualMultiArcStates.segment(
                    indexMultiArcState /*i * currentArcSize*/ + singleArcDynamicsSize_, currentArcSize - singleArcDynamicsSize_ );
            indexMultiArcState += currentArcSize;
        }

        // Reset initial multi-arc states in propagator settings and estimated parameters
        propagatorSettings_->getMultiArcPropagatorSettings( )->resetInitialStatesList( arcInitialStates );
        propagatorSettings_->updateInitialState( );
    }

    //! Function that removes the single-arc body data from propagation results before processing data
    /*!
     * Function that removes the single-arc body data from propagation results before processing data
     * \param numericalMultiArcSolution Full numerical solution of single and multi-arc bodies.
     */
    void removeSingleArcBodiesFromMultiArcSolultion( std::vector< std::map< TimeType, VectorType > >& numericalMultiArcSolution )
    {
        // Iterate over all arcs
        //        std::cout << "size numerical multi-arc solution: " << numericalMultiArcSolution.size( ) << "\n\n";
        for( unsigned int i = 0; i < numericalMultiArcSolution.size( ); i++ )
        {
            // Iterate over all times and remove single-arc bodies from solution
            for( typename std::map< TimeType, VectorType >::iterator mapIterator = numericalMultiArcSolution[ i ].begin( );
                 mapIterator != numericalMultiArcSolution[ i ].end( );
                 mapIterator++ )
            {
                VectorType fullVector = mapIterator->second;
                numericalMultiArcSolution[ i ][ mapIterator->first ] =
                        fullVector.segment( singleArcDynamicsSize_, originalMultiArcDynamicsSingleArcSize_.at( i ) );
            }
        }
    }

    //! Update original propagator settings
    void copyExtendedMultiArcInitialStatesToOriginalSettins( )
    {
        std::vector< VectorType > extendedMultiArcInitialStates =
                propagatorSettings_->getMultiArcPropagatorSettings( )->getInitialStateList( );
        std::vector< VectorType > originalMultiArcInitialStates =
                originalPopagatorSettings_->getMultiArcPropagatorSettings( )->getInitialStateList( );
        //        std::cout << "size extendedMultiArcInitialStates: " << extendedMultiArcInitialStates.size( ) << "\n\n";
        for( unsigned int i = 0; i < extendedMultiArcInitialStates.size( ); i++ )
        {
            originalMultiArcInitialStates[ i ] =
                    extendedMultiArcInitialStates.at( i ).segment( singleArcDynamicsSize_, originalMultiArcDynamicsSingleArcSize_.at( i ) );
        }

        originalPopagatorSettings_->getMultiArcPropagatorSettings( )->resetInitialStatesList( originalMultiArcInitialStates );
        originalPopagatorSettings_->updateInitialState( );
    }

    //! Object to solve multi-arc variational equations (multi-arc bodies only).
    std::shared_ptr< MultiArcVariationalEquationsSolver< StateScalarType, TimeType > > originalMultiArcSolver_;

    //! Object to solve single-arc variational equations.
    std::shared_ptr< SingleArcVariationalEquationsSolver< StateScalarType, TimeType > > singleArcSolver_;

    //! Object to solve multi-arc variational equations (single- and multi-arc bodies).
    std::shared_ptr< MultiArcVariationalEquationsSolver< StateScalarType, TimeType > > multiArcSolver_;

    //! Propagator settings, with single- and multi-arc separately
    std::shared_ptr< HybridArcPropagatorSettings< StateScalarType, TimeType > > originalPopagatorSettings_;

    //! Propagator settings, with single-arc bodies added to multi-arc list
    std::shared_ptr< HybridArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings_;

    //! Settings to be used for integrator
    std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > singleArcIntegratorSettings_;

    std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > multiArcIntegratorSettings_;

    //! Size of estimated single-arc dynamical parameters
    int singleArcDynamicsSize_;

    //! Vector containing, for each arc, the size of original estimated multi-arc dynamical parameters
    std::vector< int > originalMultiArcDynamicsSingleArcSize_;

    //! Total size of original estimated multi-arc dynamical parameters
    int originalMultiArcDynamicsSize_;

    //! Size of single arc of extended estimated multi-arc dynamical parameters
    int multiArcDynamicsSize_;

    //! Vector containing, for each arc, the total size of extended estimated multi-arc dynamical parameters
    std::vector< int > multiArcDynamicsSingleArcSize_;

    //! Estimated parameter set with single-arc dynamical parameters only
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > singleArcParametersToEstimate_;

    //! Estimated parameter set with original multi-arc dynamical parameters only
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > originalMultiArcParametersToEstimate_;

    //! Estimated parameter set with extended multi-arc dynamical parameters only
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > multiArcParametersToEstimate_;

    std::shared_ptr< HybridArcResults > variationalPropagationResults_;

    //! Times at which arcs for multi-arc solution start
    std::vector< double > arcStartTimes_;

    //! Function that retrieves the single-arc bodies' initial states as a function of time
    /*!
     *  Function that retrieves the single-arc bodies' initial states as a function of time, is used to update the multi-arc
     *  initial states after the single-arc propagation
     */
    std::function< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >( const double ) > initialStatesFromSingleArcPropagation_;

    double singleArcInitialTime_;
};

extern template class HybridArcVariationalEquationsSolver< double, double >;

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_HYBRIDARCVARIATIONALEQUATIONSSOLVER_H
