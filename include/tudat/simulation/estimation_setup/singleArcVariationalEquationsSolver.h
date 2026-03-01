/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SINGLEARCVARIATIONALEQUATIONSSOLVER_H
#define TUDAT_SINGLEARCVARIATIONALEQUATIONSSOLVER_H

#include "tudat/simulation/estimation_setup/variationalEquationsSolverBase.h"
#include "tudat/simulation/propagation_setup/singleArcDynamicsSimulator.h"

namespace tudat
{

namespace propagators
{
//! Class to manage and execute the numerical integration of variational equations of a dynamical system in a single arc.
/*!
 *  Class to manage and execute the numerical integration of variational equations of a dynamical system, in addition
 *  to the dynamics itself, in a single arc: i.e. the governing equations a single initial time, and are propagated once
 *  for the full prescribed time interval. This is in contrast to multi-arc dynamics, where the time interval is cut into
 *  pieces. In this class, the governing equations are set once, but can be re-integrated for
 *  different initial conditions using the same instance of the class.
 */
template< typename StateScalarType = double, typename TimeType = double >
class SingleArcVariationalEquationsSolver : public VariationalEquationsSolver< StateScalarType, TimeType >
{
public:
    //! Local typedefs for vector and matrix of given scalar type
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > MatrixType;
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > VectorType;

    //! Base class using statements
    using VariationalEquationsSolver< StateScalarType, TimeType >::parametersToEstimate_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::bodies_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::stateTransitionMatrixSize_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::parameterVectorSize_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::stateTransitionInterface_;

    //! Constructor
    /*!
     *  Constructor, sets up object for automatic evaluation and numerical integration of variational equations and
     *  equations of motion.
     *  \param bodies Map of bodies (with names) of all bodies in integration.
     *  \param integratorSettings Settings for numerical integrator of combined propagation of variational equations
     *  and equations of motion.
     *  \param propagatorSettings Settings for propagation of equations of motion.
     *  \param parametersToEstimate Object containing all parameters that are to be estimated and their current
     *  settings and values.
     *  \param integrateDynamicalAndVariationalEquationsConcurrently Boolean defining whether variational and dynamical
     *  equations are to be propagated concurrently (if true) or sequentially (of false)
     *  \param variationalOnlyIntegratorSettings Settings for numerical integrator when integrating only variational
     *  equations.
     *  \param clearNumericalSolution Boolean to determine whether to clear the raw numerical solution member variables
     *  (default true) after propagation and resetting of state transition interface.
     *  \param integrateEquationsOnCreation Boolean to denote whether equations should be integrated immediately at the
     *  end of this contructor (default true).
     *  \param setIntegratedResult Boolean to determine whether to automatically use the integrated results to set
     *  ephemerides (default true).
     */
    SingleArcVariationalEquationsSolver(
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< SingleArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
            const bool integrateDynamicalAndVariationalEquationsConcurrently = true,
            const bool integrateEquationsOnCreation = true ):
        VariationalEquationsSolver< StateScalarType, TimeType >(
                bodies,
                parametersToEstimate,
                propagatorSettings != nullptr ? propagatorSettings->getOutputSettingsWithCheck( )->getClearNumericalSolutions( ) : false ),
        propagatorSettings_( std::dynamic_pointer_cast< SingleArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ) )
    {
        // Check input consistency
        if( std::dynamic_pointer_cast< SingleArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ) == nullptr )
        {
            throw std::runtime_error( "Error in variational equations solver, input must be single-arc." );
        }
        else if( !checkPropagatorSettingsAndParameterEstimationConsistency< StateScalarType, TimeType >( propagatorSettings_,
                                                                                                         parametersToEstimate ) )
        {
            throw std::runtime_error(
                    "Error when making single arc variational equations solver, estimated and propagated bodies are inconsistent." );
        }

        // Create state derivative models
        std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > stateDerivativeModels =
                createStateDerivativeModels( propagatorSettings_, bodies, propagatorSettings_->getInitialTime( ) );

        // Create state derivative partials
        std::map< IntegratedStateType, orbit_determination::StateDerivativePartialsMap > stateDerivativePartials =
                simulation_setup::createStateDerivativePartials< StateScalarType, TimeType >(
                        getStateDerivativeModelMapFromVector( stateDerivativeModels ), bodies, parametersToEstimate );

        // Create object that propagates the dynamics
        dynamicsSimulator_ = std::make_shared< SingleArcDynamicsSimulator< StateScalarType, TimeType > >(
                bodies,
                propagatorSettings_,
                false,
                PredefinedSingleArcStateDerivativeModels< StateScalarType, TimeType >( stateDerivativeModels, stateDerivativePartials ) );

        // Create variational equations evaluation objects.
        variationalEquationsObject_ =
                std::make_shared< VariationalEquations >( stateDerivativePartials,
                                                          parametersToEstimate_,
                                                          dynamicsSimulator_->getDynamicsStateDerivative( )->getStateTypeStartIndices( ) );
        dynamicsSimulator_->getDynamicsStateDerivative( )->addVariationalEquations( variationalEquationsObject_ );

        // Create object that will contain and process the propagation results
        variationalPropagationResults_ = std::make_shared< SingleArcVariationalSimulationResults< StateScalarType, TimeType > >(
                dynamicsSimulator_->getSingleArcPropagationResults( ),
                this->stateTransitionMatrixSize_,
                this->parameterVectorSize_ - this->stateTransitionMatrixSize_ );

        // Integrate variational equations from initial state estimate.
        if( integrateEquationsOnCreation )
        {
            if( integrateDynamicalAndVariationalEquationsConcurrently )
            {
                integrateVariationalAndDynamicalEquations( propagatorSettings_->getInitialStates( ), true );
            }
            else
            {
                integrateVariationalAndDynamicalEquations( propagatorSettings_->getInitialStates( ), false );
            }
        }
        else
        {
            stateTransitionInterface_ = std::make_shared< SingleArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                    nullptr,
                    nullptr,
                    propagatorSettings_->getConventionalStateSize( ),
                    parameterVectorSize_,
                    variationalEquationsObject_->getStatePartialAdditionIndices( ) );
        }
    }

    SingleArcVariationalEquationsSolver(
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
            const bool integrateDynamicalAndVariationalEquationsConcurrently = true,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings =
                    std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ),
            const bool clearNumericalSolution = true,
            const bool integrateEquationsOnCreation = true,
            const bool setIntegratedResult = true,
            const bool printDependentVariableData = true,
            const bool setDependentVariablesInterface = false ):
        SingleArcVariationalEquationsSolver( bodies,
                                             validateDeprecatedSingleArcSettings( integratorSettings,
                                                                                  propagatorSettings,
                                                                                  clearNumericalSolution,
                                                                                  setIntegratedResult,
                                                                                  false,
                                                                                  printDependentVariableData,
                                                                                  false,
                                                                                  setDependentVariablesInterface ),
                                             parametersToEstimate,
                                             integrateDynamicalAndVariationalEquationsConcurrently,
                                             integrateEquationsOnCreation )
    {}

    //! Destructor
    ~SingleArcVariationalEquationsSolver( ) {}

    //! Function to integrate equations of motion only.
    /*!
     *  Function to integrate equations of motion only (in single arc).  If dynamical
     *  solution is to be processed, the environment is also updated to the new solution.
     *  \param initialStateEstimate Initial state of the equations of motion that is to be used (in same order as in
     *  parametersToEstimate_).
     */
    void integrateDynamicalEquationsOfMotionOnly( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialStateEstimate )
    {
        dynamicsSimulator_->integrateEquationsOfMotion( initialStateEstimate );
    }

    //! Function to integrate variational equations and equations of motion.
    /*!
     *  Function to integrate variational equations and equations of motion (in single arc). At the end of this function,
     *  the stateTransitionInterface_ is reset with the new state transition and sensitivity matrices. If dynamical
     *  solution is to be processed, the environment is also updated to the new solution.
     *  \param initialStateEstimate Initial state of the equations of motion that is to be used (in same order as in
     *  parametersToEstimate_).
     *  \param integrateEquationsConcurrently Variable determining whether the equations of motion are to be
     *  propagated concurrently with variational equations of motion (if true), or before variational equations (if false).
     */
    void integrateVariationalAndDynamicalEquations( const VectorType& initialStateEstimate, const bool integrateEquationsConcurrently )
    {
        if( integrateEquationsConcurrently )
        {
            // Create initial conditions from new estimate.
            MatrixType initialVariationalState =
                    this->createInitialConditions( dynamicsSimulator_->getDynamicsStateDerivative( )->convertFromOutputSolution(
                            initialStateEstimate, propagatorSettings_->getInitialTime( ) ) );

            // Propagate dynamics and variational equations
            dynamicsSimulator_->integrateEquationsOfMotion( initialVariationalState, variationalPropagationResults_ );
        }

        // Reset solution for state transition and sensitivity matrices.
        if( propagatorSettings_->getOutputSettings( )->getSetIntegratedVariationalResult( ) )
        {
            resetVariationalEquationsInterpolators( );
        }
    }

    //! Function to return the numerical solution history of numerically integrated variational equations.
    /*!
     *  Function to return the numerical solution history of numerically integrated variational equations.
     *  \return Vector of mapa of state transition matrix history (first vector entry)
     *  and sensitivity matrix history (second vector entry)
     */
    std::vector< std::map< double, Eigen::MatrixXd > > getNumericalVariationalEquationsSolution( )
    {
        std::cerr << "Warning, use of deprecated single-arc getNumericalVariationalEquationsSolution is not recommended" << std::endl;
        return std::vector< std::map< double, Eigen::MatrixXd > >(
                { getStateTransitionMatrixSolution( ), getSensitivityMatrixSolution( ) } );
    }

    std::map< double, Eigen::MatrixXd >& getStateTransitionMatrixSolution( )
    {
        return variationalPropagationResults_->getStateTransitionSolution( );
    }

    std::map< double, Eigen::MatrixXd >& getSensitivityMatrixSolution( )
    {
        return variationalPropagationResults_->getSensitivitySolution( );
    }

    const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& getEquationsOfMotionSolution( )
    {
        return dynamicsSimulator_->getEquationsOfMotionNumericalSolution( );
    }

    std::map< double, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > getEquationsOfMotionSolutionDouble( )
    {
        return dynamicsSimulator_->getEquationsOfMotionNumericalSolutionDouble( );
    }

    //! Function to return object used for numerically propagating and managing the solution of the equations of motion.
    /*!
     * Function to return object used for numerically propagating and managing the solution of the equations of motion.
     * \return Object used for numerically propagating and managing the solution of the equations of motion.
     */
    std::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > getDynamicsSimulator( )
    {
        return dynamicsSimulator_;
    }

    std::shared_ptr< VariationalEquations > getVariationalEquationsObject( )
    {
        return variationalEquationsObject_;
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

    std::shared_ptr< SingleArcVariationalSimulationResults< StateScalarType, TimeType > > getSingleArcVariationalPropagationResults( )
    {
        return variationalPropagationResults_;
    }

    std::shared_ptr< SimulationResults< StateScalarType, TimeType > > getVariationalPropagationResults( )
    {
        return getSingleArcVariationalPropagationResults( );
    }

protected:
private:
    //! Reset solutions of variational equations.
    /*!
     *  Reset solutions of variational equations (stateTransitionMatrixInterpolator_ and sensitivityMatrixInterpolator_),
     *  i.e. use numerical integration results to create new look-up tables
     *  and interpolators of state transition and sensitivity matrix through the createInterpolatorsForVariationalSolution
     *  function
     */
    void resetVariationalEquationsInterpolators( )
    {
        using namespace interpolators;
        using namespace utilities;

        // Create interpolators.
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > stateTransitionMatrixInterpolator;
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > sensitivityMatrixInterpolator;

        try
        {
            createStateTransitionAndSensitivityMatrixInterpolator( stateTransitionMatrixInterpolator,
                                                                   sensitivityMatrixInterpolator,
                                                                   variationalPropagationResults_->getStateTransitionSolution( ),
                                                                   variationalPropagationResults_->getSensitivitySolution( ),
                                                                   this->clearNumericalSolution_ );
        }
        catch( const std::exception& caughtException )
        {
            std::cerr << "Error occured when post-processing single-arc variational equation integration results, and creating "
                         "interpolators, caught error is: "
                      << std::endl
                      << std::endl;
            std::cerr << caughtException.what( ) << std::endl << std::endl;
            std::cerr << "The problem may be that there is an insufficient number of data points (epochs) at which propagation results are "
                         "produced for one or more arcs. Integrated results are given at" +
                            std::to_string( variationalPropagationResults_->getStateTransitionSolution( ).size( ) ) + " epochs"
                      << std::endl;
        }

        // Create (if non-existent) or reset state transition matrix interface
        if( stateTransitionInterface_ == nullptr )
        {
            stateTransitionInterface_ = std::make_shared< SingleArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                    stateTransitionMatrixInterpolator,
                    sensitivityMatrixInterpolator,
                    propagatorSettings_->getConventionalStateSize( ),
                    parameterVectorSize_,
                    variationalEquationsObject_->getStatePartialAdditionIndices( ) );
        }
        else
        {
            std::dynamic_pointer_cast< SingleArcCombinedStateTransitionAndSensitivityMatrixInterface >( stateTransitionInterface_ )
                    ->updateMatrixInterpolators( stateTransitionMatrixInterpolator,
                                                 sensitivityMatrixInterpolator,
                                                 variationalEquationsObject_->getStatePartialAdditionIndices( ) );
        }
    }

    //! Object used for numerically propagating and managing the solution of the equations of motion.
    std::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > dynamicsSimulator_;

    //!  Object that is used to evaluate the variational equations at the given state and time.
    std::shared_ptr< VariationalEquations > variationalEquationsObject_;

    //! Settings for propagation of equations of motion.
    std::shared_ptr< SingleArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings_;

    std::shared_ptr< SingleArcVariationalSimulationResults< StateScalarType, TimeType > > variationalPropagationResults_;
};

#if TUDAT_BUILD_EXPLICIT_INSTANTIATIONS
extern template class SingleArcVariationalEquationsSolver< double, double >;
#endif

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_SINGLEARCVARIATIONALEQUATIONSSOLVER_H
