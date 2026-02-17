/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_MULTIARCVARIATIONALEQUATIONSSOLVER_H
#define TUDAT_MULTIARCVARIATIONALEQUATIONSSOLVER_H

#include "tudat/simulation/estimation_setup/variationalEquationsSolverBase.h"

namespace tudat
{

namespace propagators
{
//! Function to transfer the initial multi-arc states from propagator settings to associated initial state estimation parameters.
/*!
 *  Function to transfer the initial multi-arc states from propagator settings to associated initial state estimation parameters.
 *  \param parametersToEstimate Full set of estimated parameters to which the initial states are to be transferred
 *  \param propagatorSettings Multi-arc propagator settings from which the initial states are to be taken
 */
template< typename StateScalarType, typename TimeType >
void setPropagatorSettingsMultiArcStatesInEstimatedDynamicalParameters(
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
        const std::shared_ptr< MultiArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings )
{
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > StateType;
    typedef std::map< std::string, std::shared_ptr< estimatable_parameters::EstimatableParameter< StateType > > > ArcWiseParameterList;

    // Get list of estimated bodies
    ArcWiseParameterList estimatedBodies =
            estimatable_parameters::getListOfBodiesWithTranslationalMultiArcStateToEstimate( parametersToEstimate );
    std::vector< std::string > bodiesWithPropagatedTranslationalState = utilities::createVectorFromMapKeys( estimatedBodies );

    std::map< std::string, unsigned int > counterArcPerBody;
    std::map< std::string, StateType > arcInitialTranslationalStates;

    // Iterate over each arc and set initial state.
    std::map< std::string, std::vector< StateType > > arcInitialTranslationalStatesVector;
    for( int arc = 0; arc < propagatorSettings->getNmberOfArcs( ); arc++ )
    {
        // Check type of dynamics
        switch( propagatorSettings->getSingleArcSettings( ).at( arc )->getStateType( ) )
        {
            case translational_state: {
                std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType, TimeType > > translationalPropagatorSettings =
                        std::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType, TimeType > >(
                                propagatorSettings->getSingleArcSettings( ).at( arc ) );

                std::vector< std::string > bodiesToIntegrate = translationalPropagatorSettings->bodiesToIntegrate_;

                // Iterate over bodies and set initial state
                for( unsigned int i = 0; i < bodiesToIntegrate.size( ); i++ )
                {
                    if( counterArcPerBody.count( bodiesToIntegrate.at( i ) ) == 0 )
                    {
                        counterArcPerBody[ bodiesToIntegrate.at( i ) ] = 0;
                        arcInitialTranslationalStatesVector[ bodiesToIntegrate.at( i ) ] = { };
                    }
                    else
                    {
                        counterArcPerBody.at( bodiesToIntegrate.at( i ) ) += 1;
                    }
                    arcInitialTranslationalStatesVector.at( bodiesToIntegrate.at( i ) )
                            .push_back( translationalPropagatorSettings->getInitialStates( ).segment( i * 6, 6 ) );
                }
                break;
            }
            default:
                std::string errorMessage = "Error, cannot yet make parameters and multi-arc propagator settings consistent for " +
                        std::to_string( propagatorSettings->getSingleArcSettings( ).at( arc )->getStateType( ) );
                throw std::runtime_error( errorMessage );
        }
    }

    for( auto itr : arcInitialTranslationalStatesVector )
    {
        arcInitialTranslationalStates[ itr.first ] = StateType( 6 * itr.second.size( ) );
        for( unsigned int k = 0; k < itr.second.size( ); k++ )
        {
            arcInitialTranslationalStates.at( itr.first ).segment( k * 6, 6 ) = itr.second.at( k );
        }
    }

    // Set information in estimation objects
    for( unsigned int i = 0; i < bodiesWithPropagatedTranslationalState.size( ); i++ )
    {
        estimatedBodies.at( bodiesWithPropagatedTranslationalState.at( i ) )
                ->setParameterValue( arcInitialTranslationalStates.at( bodiesWithPropagatedTranslationalState.at( i ) ) );
    }
}

//! Class to manage and execute the numerical integration of variational equations of a dynamical system in multiple arcs.
/*!
 *  Class to manage and execute the numerical integration of variational equations of a dynamical system, in addition
 *  to the dynamics itself,  in multiple arcs: i.e. the governing equations are propagated for a set of predescribed intervals.
 *  In this class, the governing equations are set once, but can be re-integrated for different initial conditions using the
 *  same instance of the class.
 */
template< typename StateScalarType = double, typename TimeType = double >
class MultiArcVariationalEquationsSolver : public VariationalEquationsSolver< StateScalarType, TimeType >
{
public:
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > MatrixType;
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > VectorType;
    typedef MultiArcSimulationResults< SingleArcVariationalSimulationResults, StateScalarType, TimeType > MultiArcVariationalResults;

    using VariationalEquationsSolver< StateScalarType, TimeType >::parametersToEstimate_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::bodies_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::stateTransitionMatrixSize_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::parameterVectorSize_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::stateTransitionInterface_;

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
     *  \param variationalOnlyIntegratorSettings Settings for numerical integrator when integrating only variational
     *  equations.
     *  \param clearNumericalSolution Boolean to determine whether to clear the raw numerical solution member variables
     *  (default true) after propagation and resetting of state transition interface.
     *  \param integrateEquationsOnCreation Boolean to denote whether equations should be integrated immediately at the
     *  end of this contructor (default false).
     *  \param resetMultiArcDynamicsAfterPropagation Boolean denoting whether to reset the multi-arc dynamics after
     *  propagation (default true).
     */

    MultiArcVariationalEquationsSolver(
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< MultiArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
            const bool integrateEquationsOnCreation = false ):
        VariationalEquationsSolver< StateScalarType, TimeType >(
                bodies,
                parametersToEstimate,
                propagatorSettings != nullptr ? propagatorSettings->getOutputSettingsWithCheck( )->getClearNumericalSolutions( ) : false ),
        propagatorSettings_( propagatorSettings ),
        resetMultiArcDynamicsAfterPropagation_(
                propagatorSettings != nullptr ? propagatorSettings->getOutputSettingsWithCheck( )->getSetIntegratedResult( ) : false )
    {
        if( std::dynamic_pointer_cast< MultiArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ) == nullptr )
        {
            throw std::runtime_error( "Error when making multi-arc variational equations solver, input is single-arc" );
        }
        checkMultiArcPropagatorSettingsAndParameterEstimationConsistency( propagatorSettings_,
                                                                          parametersToEstimate,
                                                                          propagatorSettings->getArcStartTimes( ),
                                                                          estimatedBodiesPerArc_,
                                                                          arcIndicesPerBody_,
                                                                          areEstimatedBodiesDifferentPerArc_ );

        if( areEstimatedBodiesDifferentPerArc_ && resetMultiArcDynamicsAfterPropagation_ )
        {
            //            throw std::runtime_error( "Error in multi-arc variational equations solver, boolean
            //            resetMultiArcDynamicsAfterPropagation should be "
            //                                      "set to false when the bodies to be estimated differ from one arc to another." );
        }

        arcWiseParametersToEstimate_.clear( );
        estimatable_parameters::getParametersToEstimatePerArcTest( parametersToEstimate,
                                                                   arcWiseParametersToEstimate_,
                                                                   propagatorSettings->getArcStartTimes( ),
                                                                   estimatedBodiesPerArc_,
                                                                   arcIndicesPerBody_ );

        parameterVectorSize_ = 0;
        stateTransitionMatrixSize_ = 0;

        for( unsigned int arc = 0; arc < estimatedBodiesPerArc_.size( ); arc++ )
        {
            arcWiseStateTransitionMatrixSize_.push_back(
                    estimatable_parameters::getSingleArcInitialDynamicalStateParameterSetSize( parametersToEstimate, arc ) );
            arcWiseParameterVectorSize_.push_back( estimatable_parameters::getSingleArcParameterSetSize( parametersToEstimate, arc ) );
        }

        dynamicsSimulator_ =
                std::make_shared< MultiArcDynamicsSimulator< StateScalarType, TimeType > >( bodies, propagatorSettings, false );

        std::vector< std::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > > singleArcDynamicsSimulators =
                dynamicsSimulator_->getSingleArcDynamicsSimulators( );

        std::vector< std::shared_ptr< SingleArcVariationalSimulationResults< StateScalarType, TimeType > > > singleArcVariationalResults;
        for( unsigned int i = 0; i < dynamicsSimulator_->getSingleArcDynamicsSimulators( ).size( ); i++ )
        {
            singleArcVariationalResults.push_back( std::make_shared< SingleArcVariationalSimulationResults< StateScalarType, TimeType > >(
                    dynamicsSimulator_->getSingleArcDynamicsSimulators( ).at( i )->getSingleArcPropagationResults( ),
                    arcWiseStateTransitionMatrixSize_.at( i ),
                    arcWiseParameterVectorSize_.at( i ) - arcWiseStateTransitionMatrixSize_.at( i ) ) );
        }
        variationalPropagationResults_ =
                std::make_shared< MultiArcSimulationResults< SingleArcVariationalSimulationResults, StateScalarType, TimeType > >(
                        singleArcVariationalResults,
                        dynamicsSimulator_->getMultiArcPropagationResults( )->getMultiArcDependentVariablesInterface( ) );

        for( unsigned int i = 0; i < singleArcDynamicsSimulators.size( ); i++ )
        {
            dynamicsStateDerivatives_.push_back( singleArcDynamicsSimulators.at( i )->getDynamicsStateDerivative( ) );

            // Create variational equations objects.
            std::map< IntegratedStateType, orbit_determination::StateDerivativePartialsMap > stateDerivativePartials =
                    simulation_setup::createStateDerivativePartials< StateScalarType, TimeType >(
                            dynamicsStateDerivatives_.at( i )->getStateDerivativeModels( ), bodies, arcWiseParametersToEstimate_[ i ] );

            std::shared_ptr< VariationalEquations > variationalEquationsObject_ =
                    std::make_shared< VariationalEquations >( stateDerivativePartials,
                                                              arcWiseParametersToEstimate_[ i ],
                                                              dynamicsStateDerivatives_.at( i )->getStateTypeStartIndices( ),
                                                              i,
                                                              arcIndicesPerBody_[ i ] );

            dynamicsStateDerivatives_.at( i )->addVariationalEquations( variationalEquationsObject_ );
        }

        numberOfArcs_ = dynamicsStateDerivatives_.size( );
        // Integrate variational equations from initial state estimate.
        if( integrateEquationsOnCreation )
        {
            integrateVariationalAndDynamicalEquations( propagatorSettings_->getInitialStateList( ), 1 );
        }
    }

    MultiArcVariationalEquationsSolver(
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
            const std::vector< double > propagationStartTimes,
            const bool integrateDynamicalAndVariationalEquationsConcurrently = true,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings =
                    std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ),
            const bool clearNumericalSolution = true,
            const bool integrateEquationsOnCreation = false,
            const bool resetMultiArcDynamicsAfterPropagation = true,
            const bool setDependentVariablesInterface = false ):
        MultiArcVariationalEquationsSolver( bodies,
                                            validateDeprecatedMultiArcSettings( integratorSettings,
                                                                                propagatorSettings,
                                                                                propagationStartTimes,
                                                                                clearNumericalSolution,
                                                                                resetMultiArcDynamicsAfterPropagation,
                                                                                setDependentVariablesInterface ),
                                            parametersToEstimate,
                                            integrateEquationsOnCreation )
    {}

    //! Destructor
    /*!
     *  Destructor
     */
    ~MultiArcVariationalEquationsSolver( ) {}

    //! Function to integrate equations of motion only.
    /*!
     *  Function to integrate equations of motion only (for all arcs).  If dynamical
     *  solution is to be processed, the environment is also updated to the new solution.
     *  \param initialStateEstimate Initial state of the equations of motion that is to be used (in same order as in
     *  parametersToEstimate_), concatenated for all arcs.
     */
    void integrateDynamicalEquationsOfMotionOnly( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialStateEstimate )
    {
        dynamicsSimulator_->integrateEquationsOfMotion( initialStateEstimate );
    }

    //! Function to integrate equations of motion only.
    /*!
     *  Function to integrate equations of motion only (for all arcs).  If dynamical
     *  solution is to be processed, the environment is also updated to the new solution.
     *  \param initialStateEstimate Initial state of the equations of motion that is to be used, as a list with separate entries
     *  for each arc.
     */
    void integrateDynamicalEquationsOfMotionOnly(
            const std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& initialStateEstimate )
    {
        dynamicsSimulator_->integrateEquationsOfMotion( initialStateEstimate );
    }

    //! Function to integrate variational equations and equations of motion.
    /*!
     *  Function to integrate variational equations and equations of motion, for all arcs. At the end of this function,
     *  the stateTransitionInterface_ is reset with the new state transition and sensitivity matrices. If dynamical
     *  solution is to be processed, the environment is also updated to the new solution.
     *  \param concatenatedInitialStates Initial state of the equations of motion that is to be used (in same order as in
     *  parametersToEstimate_), concatenated for all arcs.
     *  \param integrateEquationsConcurrently Variable determining whether the equations of motion are to be
     *  propagated concurrently with variational equations of motion (if true), or before variational equations (if false).
     */
    void integrateVariationalAndDynamicalEquations( const VectorType& concatenatedInitialStates, const bool integrateEquationsConcurrently )
    {
        std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > splitInitialState;

        int currentIndex = 0;
        for( unsigned int i = 0; i < dynamicsSimulator_->getSingleArcDynamicsSimulators( ).size( ); i++ )
        {
            int currentSize =
                    dynamicsSimulator_->getSingleArcDynamicsSimulators( ).at( i )->getPropagatorSettings( )->getConventionalStateSize( );
            splitInitialState.push_back( concatenatedInitialStates.block( currentIndex, 0, currentSize, 1 ) );
            currentIndex += currentSize;
        }

        if( currentIndex != concatenatedInitialStates.rows( ) )
        {
            throw std::runtime_error(
                    "Error when doing multi-arc variational equation integration, "
                    "input state vector size is incompatible with settings." );
        }
        integrateVariationalAndDynamicalEquations( splitInitialState, integrateEquationsConcurrently );
    }

    std::shared_ptr< MultiArcInitialStateProvider< StateScalarType > > getInitialStateProvider(
            const std::vector< VectorType >& initialStateEstimate )
    {
        std::vector< std::pair< int, int > > variationalEquationsSize =
                utilities::mergeVectorsIntoVectorOfPairs( arcWiseStateTransitionMatrixSize_, arcWiseParameterVectorSize_ );
        return std::make_shared< MultiArcInitialStateProvider< StateScalarType > >( initialStateEstimate, variationalEquationsSize );
    }

    //! Function to integrate variational equations and equations of motion.
    /*!
     *  Function to integrate variational equations and equations of motion, for all arcs. At the end of this function,
     *  the stateTransitionInterface_ is reset with the new state transition and sensitivity matrices. If dynamical
     *  solution is to be processed, the environment is also updated to the new solution.
     *  \param initialStateEstimate Initial state of the equations of motion that is to be used, in a list with initial
     *  state for each arc stored separately.
     *  \param integrateEquationsConcurrently Variable determining whether the equations of motion are to be
     *  propagated concurrently with variational equations of motion (if true), or before variational equations (if false).
     */
    void integrateVariationalAndDynamicalEquations( const std::vector< VectorType >& initialStateEstimate,
                                                    const bool integrateEquationsConcurrently )
    {
        // Propagate variational equations and equations of motion concurrently
        if( integrateEquationsConcurrently )
        {
            // Propagate dynamics and variational equations and store results in variationalPropagationResults_ object
            dynamicsSimulator_->template integrateEquationsOfMotion< MultiArcVariationalResults >(
                    variationalPropagationResults_, getInitialStateProvider( initialStateEstimate ) );

            // Ensure consistency between parameters and propagator settings
            setPropagatorSettingsMultiArcStatesInEstimatedDynamicalParameters< StateScalarType, TimeType >( parametersToEstimate_,
                                                                                                            propagatorSettings_ );
        }

        // Reset solution for state transition and sensitivity matrices.
        if( propagatorSettings_->getOutputSettings( )->getSetIntegratedVariationalResult( ) )
        {
            resetVariationalEquationsInterpolators( );
        }
    }

    //! Function to return object used for numerically propagating and managing the solution of the equations of motion.
    /*!
     * Function to return object used for numerically propagating and managing the solution of the equations of motion.
     * \return Object used for numerically propagating and managing the solution of the equations of motion.
     */
    std::shared_ptr< MultiArcDynamicsSimulator< StateScalarType, TimeType > > getDynamicsSimulator( )
    {
        return dynamicsSimulator_;
    }

    //! Function to retrieve the dynamics simulator object (as base-class pointer)
    /*!
     * Function to retrieve the dynamics simulator object (as base-class pointer)
     * \return Dynamics simulator object (as base-class pointer)
     */
    std::shared_ptr< DynamicsSimulator< StateScalarType, TimeType > > getDynamicsSimulatorBase( )
    {
        return getDynamicsSimulator( );
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
        simulation_setup::setInitialStateVectorFromParameterSet< StateScalarType, TimeType >( parametersToEstimate_, propagatorSettings_ );

        for( int i = 0; i < numberOfArcs_; i++ )
        {
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > newParametersValues =
                    propagatorSettings_->getSingleArcSettings( ).at( i )->getInitialStates( );
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > arcWiseParametersValues =
                    arcWiseParametersToEstimate_.at( i )->template getFullParameterValues< StateScalarType >( );
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > newArcWiseParametersValues = arcWiseParametersValues;
            newArcWiseParametersValues.segment( 0, newParametersValues.size( ) ) = newParametersValues;
            arcWiseParametersToEstimate_.at( i )->resetParameterValues( newArcWiseParametersValues );
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

    //! Function to return the numerical solution history of integrated variational equations, per arc.
    /*!
     *  Function to return the numerical solution history of integrated variational equations, per arc.
     *  Each vector entry contains the results of a single arc, stored in a vector of maps. Inner vector has size two: first entry
     *  is state transition matrix history, second is sensitivity matrix history, both stored as maps. Key of map denotes time,
     *  values are matrices.
     *  \return Numerical solution history of integrated variational equations, per arc.
     */
    std::vector< std::vector< std::map< double, Eigen::MatrixXd > > > getNumericalVariationalEquationsSolution( )
    {
        std::vector< std::vector< std::map< double, Eigen::MatrixXd > > > fullSolution;
        fullSolution.resize( numberOfArcs_ );
        for( int i = 0; i < numberOfArcs_; i++ )
        {
            fullSolution[ i ].push_back( variationalPropagationResults_->getSingleArcResults( ).at( i )->getStateTransitionSolution( ) );
            fullSolution[ i ].push_back( variationalPropagationResults_->getSingleArcResults( ).at( i )->getSensitivitySolution( ) );
        }
        std::cerr << "Warning, use of deprecated multi-arc getNumericalVariationalEquationsSolution is not recommended" << std::endl;
        return fullSolution;
    }

    //! Function to return list of start times of each arc. NOTE: This list is updated after every propagation.
    /*!
     * Function to return list of start times of each arc. NOTE: This list is updated after every propagation.
     * \return List of start times of each arc. NOTE: This list is updated after every propagation.
     */
    //    std::vector< double > getArcStartTimes( )
    //    {
    //        return arcStartTimes_;
    //    }

    std::vector< std::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > > getDynamicsStateDerivatives( )
    {
        return dynamicsStateDerivatives_;
    }

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > > getArcWiseParametersToEstimate( )
            const
    {
        return arcWiseParametersToEstimate_;
    }

    std::shared_ptr< MultiArcVariationalResults > getMultiArcVariationalPropagationResults( )
    {
        return variationalPropagationResults_;
    }

    std::shared_ptr< SimulationResults< StateScalarType, TimeType > > getVariationalPropagationResults( )
    {
        return getMultiArcVariationalPropagationResults( );
    }

protected:
private:
    //! Reset solutions of variational equations.
    /*!
     *  Reset solutions of variational equations (stateTransitionMatrixInterpolator_ and sensitivityMatrixInterpolator_) for each
     *  arc,*  i.e. use numerical integration results to create new look-up tables and interpolators of state transition and
     *  sensitivity matrices.
     */
    void resetVariationalEquationsInterpolators( )
    {
        using namespace interpolators;

        // Allocate interpolator vectors
        std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
                stateTransitionMatrixInterpolators;
        std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
                sensitivityMatrixInterpolators;
        stateTransitionMatrixInterpolators.resize( variationalPropagationResults_->getSingleArcResults( ).size( ) );
        sensitivityMatrixInterpolators.resize( variationalPropagationResults_->getSingleArcResults( ).size( ) );

        // Create interpolators.
        std::vector< double > arcStartTimesToUse;
        std::vector< double > arcEndTimesToUse;

        for( unsigned int i = 0; i < variationalPropagationResults_->getSingleArcResults( ).size( ); i++ )
        {
            arcStartTimesToUse.push_back( dynamicsSimulator_->getArcStartTimes( ).at( i ) );
            arcEndTimesToUse.push_back( dynamicsSimulator_->getArcEndTimes( ).at( i ) );

            try
            {
                createStateTransitionAndSensitivityMatrixInterpolator(
                        stateTransitionMatrixInterpolators[ i ],
                        sensitivityMatrixInterpolators[ i ],
                        variationalPropagationResults_->getSingleArcResults( ).at( i )->getStateTransitionSolution( ),
                        variationalPropagationResults_->getSingleArcResults( ).at( i )->getSensitivitySolution( ),
                        this->clearNumericalSolution_ );
            }
            catch( const std::exception& caughtException )
            {
                std::cerr << "Error occured when post-processing multi-arc variational equation integration results, and creating "
                             "interpolators in arc" +
                                std::to_string( i ) + ", caught error is: "
                          << std::endl
                          << std::endl;
                std::cerr << caughtException.what( ) << std::endl << std::endl;
                std::cerr << "The problem may be that there is an insufficient number of data points (epochs) at which propagation results "
                             "are produced for one or more arcs. Integrated results are given at" +
                                std::to_string( variationalPropagationResults_->getSingleArcResults( )
                                                        .at( 0 )
                                                        ->getStateTransitionSolution( )
                                                        .size( ) ) +
                                " epochs"
                          << std::endl;
            }
        }

        // Create stare transition matrix interface if needed, reset otherwise.
        if( stateTransitionInterface_ == nullptr )
        {
            stateTransitionInterface_ = std::make_shared< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface< StateScalarType > >(
                    stateTransitionMatrixInterpolators,
                    sensitivityMatrixInterpolators,
                    propagatorSettings_->getArcStartTimes( ),
                    arcStartTimesToUse,
                    arcEndTimesToUse,
                    parametersToEstimate_,
                    propagatorSettings_->getSingleArcSettings( ).at( 0 )->getConventionalStateSize( ),
                    parametersToEstimate_->getParameterSetSize( ),
                    getArcWiseStatePartialAdditionIndices( ) );
        }
        else
        {
            std::dynamic_pointer_cast< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface< StateScalarType > >(
                    stateTransitionInterface_ )
                    ->updateMatrixInterpolators( stateTransitionMatrixInterpolators,
                                                 sensitivityMatrixInterpolators,
                                                 arcStartTimesToUse,
                                                 arcEndTimesToUse,
                                                 getArcWiseStatePartialAdditionIndices( ) );
        }
    }

    std::vector< std::vector< std::pair< int, int > > > getArcWiseStatePartialAdditionIndices( )
    {
        std::vector< std::vector< std::pair< int, int > > > partialIndices;
        for( unsigned int i = 0; i < dynamicsStateDerivatives_.size( ); i++ )
        {
            partialIndices.push_back(
                    dynamicsStateDerivatives_.at( i )->getVariationalEquationsCalculator( )->getStatePartialAdditionIndices( ) );
        }
        return partialIndices;
    }

    //! Object to propagate the dynamics for all arcs.
    std::shared_ptr< MultiArcDynamicsSimulator< StateScalarType, TimeType > > dynamicsSimulator_;

    //    //! Numerical solution history of integrated variational equations, per arc.
    //    /*!
    //     *  Numerical solution history of integrated variational equations, per arc.
    //     *  Each vector entry contains the results of a single arc, stored in a vector of maps. Inner vector has size two: first entry
    //     *  is state transition matrix history, second is sensitivity matrix history, both stored as maps. Key of map denotes time,
    //     *  values are matrices.
    //     */
    //    std::vector< std::vector< std::map< double, Eigen::MatrixXd > > > variationalEquationsSolution_;

    //    //! List of start times of each arc. NOTE: This list is updated after every propagation.
    //    std::vector< double > arcStartTimes_;

    //    std::vector< double > arcEndTimes_;

    //! Settings for propagation of equations of motion.
    std::shared_ptr< MultiArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings_;

    //! State derivative models for each arc (retrieved from dynamicsSimulator_).
    std::vector< std::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > > dynamicsStateDerivatives_;

    //! Number of arcs over which propagation is to be performed.
    int numberOfArcs_;

    //! Boolean denoting whether to reset the multi-arc dynamics after propagation.
    bool resetMultiArcDynamicsAfterPropagation_;

    //! Map containing, for each arc, a vector with the names of the bodies whose initial states are to be estimated.
    std::map< int, std::vector< std::string > > estimatedBodiesPerArc_;

    //! Map containing, for each arc, a map where the keys are the propagated bodies and the elements give the arc index body-wise
    //! (e.g. arc j can be the kth arc for which body i is propagated).
    std::map< int, std::map< std::string, int > > arcIndicesPerBody_;

    //! Vector containing the size of the state transition matrix, for each arc.
    std::vector< int > arcWiseStateTransitionMatrixSize_;

    //! Vector containing the size of the sensitivity matrix, for each arc.
    std::vector< int > arcWiseParameterVectorSize_;

    //! Vector with arc-wise parameters to be estimated.
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > > arcWiseParametersToEstimate_;

    //! Boolean denoting whether the estimated bodies are different from one arc to another.
    bool areEstimatedBodiesDifferentPerArc_;

    std::shared_ptr< MultiArcVariationalResults > variationalPropagationResults_;
};

extern template class MultiArcVariationalEquationsSolver< double, double >;

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_MULTIARCVARIATIONALEQUATIONSSOLVER_H
