/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_VARIATIONALEQUATIONSSOLVER_H
#define TUDAT_VARIATIONALEQUATIONSSOLVER_H

#include <boost/make_shared.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>

#include "Tudat/Basics/utilities.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Mathematics/Interpolators/interpolator.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/Propagators/stateTransitionMatrixInterface.h"
#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"
#include "Tudat/Astrodynamics/Ephemerides/tabulatedEphemeris.h"
#include "Tudat/SimulationSetup/EstimationSetup/createStateDerivativePartials.h"
#include "Tudat/SimulationSetup/EstimationSetup/createEstimatableParameters.h"

namespace tudat
{

namespace propagators
{


//! Base class to manage and execute the numerical integration of equations of motion and variational equations.
/*!
 *  Base class to manage and execute the numerical integration of equations of motion and variational equations.
 *  Governing equations are set once, but can be re-integrated for different initial conditions using the same
 *  instance of the class. Derived classes define the specific kind of integration that is performed
 *  (single-arc/multi-arc; dynamics/variational equations, etc.)
 */
template< typename StateScalarType = double, typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< StateScalarType, TimeType >::value, int >::type = 0 >
class VariationalEquationsSolver
{
public:

    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > MatrixType;
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > VectorType;

    //! Constructor
    /*!
     *  Constructor, sets up object for automatic evaluation and numerical integration of variational equations and
     *  equations of motion.
     *  \param bodyMap Map of bodies (with names) of all bodies in integration.
     *  \param parametersToEstimate Object containing all parameters that are to be estimated and their current
     *  settings and values.
     *  \param clearNumericalSolution Boolean to determine whether to clear the raw numerical solution member variables
     *  (default true) after propagation and resetting of state transition interface.
     */
    VariationalEquationsSolver(
            const simulation_setup::NamedBodyMap& bodyMap,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
            const bool clearNumericalSolution = 1 ):
        parametersToEstimate_( parametersToEstimate ),
        bodyMap_( bodyMap ),
        stateTransitionMatrixSize_( parametersToEstimate_->getInitialDynamicalStateParameterSize( ) ),
        parameterVectorSize_( parametersToEstimate_->getParameterSetSize( ) ),
        clearNumericalSolution_( clearNumericalSolution )
    { }

    //! Destructor
    virtual ~VariationalEquationsSolver( ){ }

    //! Pure virtual function to integrate variational equations and equations of motion.
    /*!
     *  Pure virtual function to integrate variational equations and equations of motion, to be implemented in derived
     *  class
     *  \param initialStateEstimate Initial state of the equations of motion that is to be used.
     *  \param integrateEquationsConcurrently Variable determining whether the equations of motion are to be
     *  propagated concurrently with variational equations of motion (if true), or before variational equations (if false).
     */
    virtual void integrateVariationalAndDynamicalEquations(
            const VectorType& initialStateEstimate, const bool integrateEquationsConcurrently ) = 0;

    //! Pure virtual function to integrate equations of motion only.
    /*!
     *  Pure virtual function to integrate equations of motion only, to be implemented in derived
     *  class
     *  \param initialStateEstimate Initial state of the equations of motion that is to be used.
     */
    virtual void integrateDynamicalEquationsOfMotionOnly(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialStateEstimate ) = 0;


    //! Function to get the list of objects representing the parameters that are to be integrated.
    /*!
     *  Function to get the list of objects representing the parameters that are to be integrated.
     *  \return List of objects representing the parameters that are to be integrated.
     */
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > getParametersToEstimate( )
    {
        return parametersToEstimate_;
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
    virtual void resetParameterEstimate( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > newParameterEstimate,
                                         const bool areVariationalEquationsToBeIntegrated = true ) = 0;

    //! Function to get the state transition matric interface object.
    /*!
     *  Function to get the state transition matric interface object.
     *  \return The state transition matric interface object.
     */
    std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > getStateTransitionMatrixInterface( )
    {
        return stateTransitionInterface_;
    }

    //! Pure virtual function to retrieve the dynamics simulator object (as base-class pointer)
    /*!
     * Pure virtual function to retrieve the dynamics simulator object (as base-class pointer)
     * \return Dynamics simulator object (as base-class pointer)
     */
    virtual std::shared_ptr< DynamicsSimulator< StateScalarType, TimeType > > getDynamicsSimulatorBase( ) = 0;



protected:


    //! Create initial matrix of numerical soluation to variational + dynamical equations.
    /*!
     *  Create initial matrix of numerical soluation to variational + dynamical equations. The structure of the matrix is
     *  [Phi;S;y], with Phi the state transition matrix, S the sensitivity matrix y the state vector.
     *  \param initialStateEstimate vector of initial state (position/velocity) of bodies to be integrated numerically.
     *  order determined by order of bodiesToIntegrate_.
     *  \return Initial matrix of numerical soluation to variation + state equations.
     */
    MatrixType createInitialConditions( const VectorType initialStateEstimate )
    {
        if( stateTransitionMatrixSize_ != initialStateEstimate.rows( ) )
        {
            throw std::runtime_error( "Error when getting initial condition for variational equations, sizes are incompatible." );
        }

        // Initialize initial conditions to zeros.
        MatrixType varSystemInitialState = MatrixType( stateTransitionMatrixSize_,
                                                       parameterVectorSize_ + 1 ).setZero( );

        // Set initial state transition matrix to identity
        varSystemInitialState.block( 0, 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_ ).setIdentity( );

        // Set initial body states to current estimate of initial body states.
        varSystemInitialState.block( 0, parameterVectorSize_,
                                     stateTransitionMatrixSize_, 1 ) = initialStateEstimate;

        return varSystemInitialState;
    }

    //! Create initial matrix of numerical soluation to variational equations
    /*!
     *  Create initial matrix of numerical soluation to variational equations, with structure [Phi;S]. Initial state
     *  transition matrix Phi is identity matrix. Initial sensitivity matrix S is all zeros.
     *  \return Initial matrix solution to variational equations.
     */
    Eigen::MatrixXd createInitialVariationalEquationsSolution( )
    {
        // Initialize initial conditions to zeros.
        Eigen::MatrixXd varSystemInitialState = Eigen::MatrixXd::Zero(
                    stateTransitionMatrixSize_, parameterVectorSize_ );

        // Set initial state transition matrix to identity
        varSystemInitialState.block( 0, 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_ ).setIdentity( );

        return varSystemInitialState;
    }

    //! Object containing all parameters that are to be estimated and their current  settings and values.
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate_ ;

    //! Map of bodies (with names) of all bodies in integration.
    simulation_setup::NamedBodyMap bodyMap_;

    //! Size (rows and columns are equal) of state transition matrix.
    int stateTransitionMatrixSize_;

    //! Number of rows in sensitivity matrix
    int parameterVectorSize_;

    //! Boolean to determine whether to clear the raw numerical solution member variables after propagation
    /*!
     *  Boolean to determine whether to clear the raw numerical solution member variables after propagation
     *  and resetting of state transition interface.
     */
    bool clearNumericalSolution_;

    //! Object used for interpolating numerical results of state transition and sensitivity matrix.
    std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface_;
};

//! Function to separate the time histories of the sensitivity and state transition matrices from a full numerical solution.
/*!
 *  Function to separate the time histories of the sensitivity and state transition matrices from a full numerical solution,
 *  in which the solution is represented as a single matrix block per time value.
 *  NOTE: numericalIntegrationResult contents are deleted by this function (all information is conserved in
 *  variationalEquationsSolution.
 *  \param numericalIntegrationResult Full time history from which separate matrix histories are to be retrieved.
 *  \param variationalEquationsSolution Vector of two matrix histories (returned by reference). First vector entry
 *  is state transition matrix history, second entry is sensitivity matrix history.
 *  \param stateTransitionStartIndices First row and column (first and second) of state transition matrix in entries of
 *  numericalIntegrationResult.
 *  \param sensitivityStartIndices First row and column (first and second) of sensitivity matrix in entries of
 *  numericalIntegrationResult.
 *  \param stateTransitionMatrixSize Size (rows and columns are equal) of state transition matrix.
 *  \param parameterSetSize Number of rows in sensitivity matrix
 */
template< typename TimeType, typename StateScalarType >
void setVariationalEquationsSolution(
        std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >&
        numericalIntegrationResult,
        std::vector< std::map< double, Eigen::MatrixXd > >& variationalEquationsSolution,
        const std::pair< int, int > stateTransitionStartIndices,
        const std::pair< int, int > sensitivityStartIndices,
        const int stateTransitionMatrixSize,
        const int parameterSetSize )
{
    variationalEquationsSolution.clear( );
    variationalEquationsSolution.resize( 2 );

    for( typename std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >::iterator
         integrationIterator = numericalIntegrationResult.begin( );
         integrationIterator != numericalIntegrationResult.end( ); )
    {
        // Set result for state transition matrix in each time step.
        variationalEquationsSolution[ 0 ][ integrationIterator->first ] =
                ( integrationIterator->second.block( stateTransitionStartIndices.first, stateTransitionStartIndices.second,
                                                     stateTransitionMatrixSize,
                                                     stateTransitionMatrixSize ) ).template cast< double >( );

        // Set result for sensitivity matrix in each time step.
        variationalEquationsSolution[ 1 ][ integrationIterator->first ] =
                ( integrationIterator->second.block( sensitivityStartIndices.first, sensitivityStartIndices.second,
                                                     stateTransitionMatrixSize,
                                                     parameterSetSize -
                                                     stateTransitionMatrixSize ) ).template cast< double >( );
        numericalIntegrationResult.erase( integrationIterator++ );
    }
}

//! Function to create interpolators for state transition and sensitivity matrices from numerical results.
/*!
 * Function to create interpolators for state transition and sensitivity matrices from numerical results.
 * \param stateTransitionMatrixInterpolator Interpolator object for state transition matrix (returned by reference).
 * \param sensitivityMatrixInterpolator Interpolator object for sensitivity matrix (returned by reference).
 * \param variationalEquationsSolution Vector of two matrix histories. First vector entry
 *  is state transition matrix history, second entry is sensitivity matrix history.
 * \param clearRawSolution Boolean denoting whether to clear entries of variationalEquationsSolution after creation
 * of interpolators.
 */
void createStateTransitionAndSensitivityMatrixInterpolator(
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >&
        stateTransitionMatrixInterpolator,
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >&
        sensitivityMatrixInterpolator,
        std::vector< std::map< double, Eigen::MatrixXd > >& variationalEquationsSolution,
        const bool clearRawSolution = 1 );

//! Function to check the consistency between propagation settings of equations of motion, and estimated parameters.
/*!
 *  Function to check the consistency between propagation settings of equations of motion, and estimated parameters.
 *  In particular, it is presently required that the set of propagated states is equal to the set of estimated states.
 *  \param propagatorSettings Settings for propagation of equations of motion.
 *  \param parametersToEstimate Object containing all parameters that are to be estimated and their current
 *  settings and values.
 *  \return True if settings are consistent
 */
template< typename StateScalarType = double, typename TimeType = double >
bool checkPropagatorSettingsAndParameterEstimationConsistency(
        const std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate )
{
    bool isInputConsistent = 1;

    // Check type of dynamics
    switch( propagatorSettings->getStateType( ) )
    {
    case translational_state:
    {
        std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > translationalPropagatorSettings =
                std::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType > >( propagatorSettings );

        // Retrieve estimated and propagated translational states, and check equality.
        std::vector< std::string > propagatedBodies = translationalPropagatorSettings->bodiesToIntegrate_;
        std::vector< std::string > estimatedBodies = estimatable_parameters::getListOfBodiesWithTranslationalStateToEstimate(
                    parametersToEstimate );

        if( static_cast< unsigned int >( translationalPropagatorSettings->getPropagatedStateSize( ) ) !=
                propagatedBodies.size( ) * 6 )
        {
            throw std::runtime_error( "Error when propagating variational equations, tranbbslational state vectors not of size 6." );
        }

        if( propagatedBodies.size( ) != estimatedBodies.size( ) )
        {
            std::string errorMessage = "Error, propagated and estimated body vector sizes are inconsistent " +
                    std::to_string( propagatedBodies.size( ) ) + " " +
                    std::to_string( estimatedBodies.size( ) );
            throw std::runtime_error( errorMessage );
            isInputConsistent = 0;
        }
        else
        {
            for( unsigned int i = 0; i < propagatedBodies.size( ); i++ )
            {
                if( propagatedBodies.at( i ) != estimatedBodies.at( i ) )
                {
                    std::string errorMessage = "Error, propagated and estimated body vectors inconsistent at index" +
                            std::string( propagatedBodies.at( i ) ) + " " +
                            std::string( estimatedBodies.at( i ) );
                    throw std::runtime_error( errorMessage );
                    isInputConsistent = 0;
                }
            }

        }
        break;
    }
    case rotational_state:
    {
        std::shared_ptr< RotationalStatePropagatorSettings< StateScalarType > > rotationalPropagatorSettings =
                std::dynamic_pointer_cast< RotationalStatePropagatorSettings< StateScalarType > >( propagatorSettings );

        // Retrieve estimated and propagated translational states, and check equality.
        std::vector< std::string > propagatedBodies = rotationalPropagatorSettings->bodiesToIntegrate_;
        std::vector< std::string > estimatedBodies = estimatable_parameters::getListOfBodiesWithRotationalStateToEstimate(
                    parametersToEstimate );
        if( propagatedBodies.size( ) != estimatedBodies.size( ) )
        {
            std::string errorMessage = "Error, propagated and estimated body vector sizes are inconsistent " +
                    std::to_string( propagatedBodies.size( ) ) + " " +
                    std::to_string( estimatedBodies.size( ) );
            throw std::runtime_error( errorMessage );
            isInputConsistent = 0;
        }
        else
        {
            for( unsigned int i = 0; i < propagatedBodies.size( ); i++ )
            {
                if( propagatedBodies.at( i ) != estimatedBodies.at( i ) )
                {
                    std::string errorMessage = "Error, propagated and estimated body vectors inconsistent at index" +
                            std::string( propagatedBodies.at( i ) ) + " " +
                            std::string( estimatedBodies.at( i ) );
                    throw std::runtime_error( errorMessage );
                    isInputConsistent = 0;
                }
            }

        }
        break;
    }
    case hybrid:
    {
        std::shared_ptr< MultiTypePropagatorSettings< StateScalarType > > multiTypePropagatorSettings =
                std::dynamic_pointer_cast< MultiTypePropagatorSettings< StateScalarType > >( propagatorSettings );
        isInputConsistent = true;


        for( auto settingIterator = multiTypePropagatorSettings->propagatorSettingsMap_.begin( );
             settingIterator !=  multiTypePropagatorSettings->propagatorSettingsMap_.end( ); settingIterator++ )
        {
            for( unsigned int i = 0; i < settingIterator->second.size( ); i++ )
            {
                if( !checkPropagatorSettingsAndParameterEstimationConsistency(
                            settingIterator->second.at( i ), parametersToEstimate ) )
                {
                    isInputConsistent = false;
                }
            }
        }

        if( estimatable_parameters::getListOfBodiesWithTranslationalStateToEstimate(
                    parametersToEstimate ).size( ) > 0 && multiTypePropagatorSettings->propagatorSettingsMap_.count( translational_state ) == 0 )
        {
            throw std::runtime_error( "Error, estimating but not propagating translational dynamics" );
            isInputConsistent = false;
        }


        if( estimatable_parameters::getListOfBodiesWithRotationalStateToEstimate(
                    parametersToEstimate ).size( ) > 0 && multiTypePropagatorSettings->propagatorSettingsMap_.count( rotational_state ) == 0 )
        {
            throw std::runtime_error( "Error, estimating but not propagating rotational dynamics" );
            isInputConsistent = false;
        }
        break;
    }
    default:
        std::string errorMessage = "Error, cannot yet check consistency of propagator settings for type " +
                std::to_string( propagatorSettings->getStateType( ) );
        throw std::runtime_error( errorMessage );
    }
    return isInputConsistent;
}

//! Function to check the consistency between multi-arc propagation settings of equations of motion, and estimated parameters.
/*!
 *  Function to check the consistency between multi-arc propagation settings of equations of motion, and estimated parameters.
 *  In particular, it is presently required that the set of propagated states is equal to the set of estimated states.
 *  \param propagatorSettings Settings for propagation of equations of motion.
 *  \param parametersToEstimate Object containing all parameters that are to be estimated and their current
 *  settings and values.
 *  \param arcStartTimes Times at which the dynamics arcs start
 *  \return True if settings are consistent
 */
template< typename StateScalarType = double, typename TimeType = double >
bool checkMultiArcPropagatorSettingsAndParameterEstimationConsistency(
        const std::shared_ptr< MultiArcPropagatorSettings< StateScalarType > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
        const std::vector< double > arcStartTimes )
{
    bool isInputConsistent = 1;

    // Get list of objets and associated bodies to estimate initial arc-wise translational states
    typedef std::map< std::string, std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > > ArcWiseParameterList;
    ArcWiseParameterList estimatedBodies = estimatable_parameters::getListOfBodiesWithTranslationalMultiArcStateToEstimate(
                parametersToEstimate );

    // Iterate over all parameters and check consistency
    for( typename ArcWiseParameterList::const_iterator parameterIterator = estimatedBodies.begin( ); parameterIterator !=
         estimatedBodies.end( ); parameterIterator++ )
    {
        // Get arc start times of current parameter
        std::vector< double > parameterArcStartTimes =
                std::dynamic_pointer_cast< estimatable_parameters::
                ArcWiseInitialTranslationalStateParameter< StateScalarType > >(
                    parameterIterator->second )->getArcStartTimes( );

        // Check if arc times are (almost) exactly the same
        if( arcStartTimes.size( ) != parameterArcStartTimes.size( ) )
        {
            isInputConsistent = false;
            throw std::runtime_error( "Error, arc times for " + parameterIterator->first + " have incompatible size with estimation" );
        }
        else
        {
            for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
            {
                if( std::fabs( arcStartTimes.at( i ) - parameterArcStartTimes.at( i ) ) >
                        std::max( 4.0 * parameterArcStartTimes.at( i ) * std::numeric_limits< double >::epsilon( ), 1.0E-12 ) )
                {
                    isInputConsistent = false;
                    throw std::runtime_error( "Error, arc time for " + parameterIterator->first + " is incompatible with estimation" );
                }
            }
        }
    }

    std::map< IntegratedStateType, std::vector< std::string > > propagatedStateTypes;

    // Iterate over each arc in propagator settings and check consistency
    for( int arc = 0; arc < propagatorSettings->getNmberOfArcs( ); arc++ )
    {
        // Check type of dynamics
        switch( propagatorSettings->getSingleArcSettings( ).at( arc )->getStateType( ) )
        {
        case translational_state:
        {
            std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > translationalPropagatorSettings =
                    std::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType > >(
                        propagatorSettings->getSingleArcSettings( ).at( arc ) );

            // Retrieve estimated and propagated translational states, and check equality.
            std::vector< std::string > propagatedBodies = translationalPropagatorSettings->bodiesToIntegrate_;
            if( arc == 0 )
            {
                propagatedStateTypes[ translational_state ] = propagatedBodies;
            }
            else
            {
                if( propagatedBodies.size( ) != propagatedStateTypes.at( translational_state ).size( ) )
                {
                    isInputConsistent = false;
                    std::string errorMessage = "Error, propagated body vector sizes are inconsistent between arcs " +
                            std::to_string( propagatedBodies.size( ) ) + " " +
                            std::to_string( propagatedStateTypes[ translational_state ].size( ) ) +
                            " when checking multi-arc estimation/propagation consistency";
                    throw std::runtime_error( errorMessage );
                }
                else
                {
                    for( unsigned int i = 0; i < propagatedBodies.size( ); i++ )
                    {
                        if( propagatedBodies.at( i ) != propagatedStateTypes[ translational_state ].at( i ) )
                        {
                            isInputConsistent = false;
                            std::string errorMessage = "Error, propagated body vector sizes are inconsistent between arcs at index  " +
                                    std::to_string( i ) + " " +
                                    std::string( propagatedBodies.at( i ) ) + " " +
                                    std::string( propagatedStateTypes[ translational_state ].at( i ) ) +
                                    " when checking multi-arc estimation/propagation consistency";
                            throw std::runtime_error( errorMessage );
                        }
                    }

                }
            }
            break;
        }
        default:
            std::string errorMessage = "Error, cannot yet check consistency of multi-arc propagator settings for type " +
                    std::to_string( propagatorSettings->getSingleArcSettings( ).at( arc )->getStateType( ) );
            throw std::runtime_error( errorMessage );
        }
    }

    if( estimatedBodies.size( ) != propagatedStateTypes[ translational_state ].size( ) )
    {
        isInputConsistent = false;
        std::string errorMessage = "Error, propagated body vector sizes are inconsistent " +
                std::to_string( propagatedStateTypes[ translational_state ].size( ) ) + " " +
                std::to_string( estimatedBodies.size( ) ) +
                " when checking multi-arc estimation/propagation consistency";
        throw std::runtime_error( errorMessage );

        for( unsigned int i = 0; i < propagatedStateTypes[ translational_state ].size( ); i++ )
        {
            if( estimatedBodies.count( propagatedStateTypes[ translational_state ].at( i ) ) == 0 )
            {
                isInputConsistent = false;
                std::string errorMessage = "Error, propagated body " +
                        std::string( propagatedStateTypes[ translational_state ].at( i ) ) + " " +
                        " not found in estimated body list when checking multi-arc estimation/propagation consistency";
                throw std::runtime_error( errorMessage );
            }
        }
    }

    return isInputConsistent;
}

//! Class to manage and execute the numerical integration of variational equations of a dynamical system in a single arc.
/*!
 *  Class to manage and execute the numerical integration of variational equations of a dynamical system, in addition
 *  to the dynamics itself, in a single arc: i.e. the governing equations a single initial time, and are propagated once
 *  for the full prescribed time interval. This is in contrast to multi-arc dynamics, where the time interval is cut into
 *  pieces. In this class, the governing equations are set once, but can be re-integrated for
 *  different initial conditions using the same instance of the class.
 */
template< typename StateScalarType = double, typename TimeType = double >
class SingleArcVariationalEquationsSolver: public VariationalEquationsSolver< StateScalarType, TimeType >
{
public:

    //! Local typedefs for vector and matrix of given scalar type
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > MatrixType;
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > VectorType;

    //! Base class using statements
    using VariationalEquationsSolver< StateScalarType, TimeType >::parametersToEstimate_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::bodyMap_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::stateTransitionMatrixSize_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::parameterVectorSize_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::stateTransitionInterface_;

    //! Constructor
    /*!
     *  Constructor, sets up object for automatic evaluation and numerical integration of variational equations and
     *  equations of motion.
     *  \param bodyMap Map of bodies (with names) of all bodies in integration.
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
            const simulation_setup::NamedBodyMap& bodyMap,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
            const bool integrateDynamicalAndVariationalEquationsConcurrently = true,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings
            = std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ),
            const bool clearNumericalSolution = true,
            const bool integrateEquationsOnCreation = true,
            const bool setIntegratedResult = true ):
        VariationalEquationsSolver< StateScalarType, TimeType >(
            bodyMap, parametersToEstimate, clearNumericalSolution ),
        integratorSettings_( integratorSettings ),
        propagatorSettings_( std::dynamic_pointer_cast< SingleArcPropagatorSettings< StateScalarType > >(propagatorSettings ) ),
        variationalOnlyIntegratorSettings_( variationalOnlyIntegratorSettings )
    {
        if( std::dynamic_pointer_cast< SingleArcPropagatorSettings< StateScalarType >  >( propagatorSettings ) == nullptr )
        {
            throw std::runtime_error( "Error in variational equations solver, input must be single-arc." );
        }

        // Check input consistency
        if( !checkPropagatorSettingsAndParameterEstimationConsistency< StateScalarType, TimeType >(
                    propagatorSettings_, parametersToEstimate ) )
        {
            throw std::runtime_error(
                        "Error when making single arc variational equations solver, estimated and propagated bodies are inconsistent." );
        }
        else
        {

            std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > stateDerivativeModels =
                    createStateDerivativeModels( propagatorSettings_, bodyMap, integratorSettings_->initialTime_ );

            // Create state derivative partials
            std::map< IntegratedStateType, orbit_determination::StateDerivativePartialsMap >
                    stateDerivativePartials =
                    simulation_setup::createStateDerivativePartials
                    < StateScalarType, TimeType >(
                        getStateDerivativeModelMapFromVector( stateDerivativeModels ), bodyMap, parametersToEstimate );

            // Create simulation object for dynamics only.
            if( propagatorSettings_->getDependentVariablesToSave( ) != nullptr )
            {
                propagatorSettings_->getDependentVariablesToSave( )->stateDerivativePartials_ = stateDerivativePartials;
            }

            dynamicsSimulator_ = std::make_shared< SingleArcDynamicsSimulator< StateScalarType, TimeType > >(
                        bodyMap, integratorSettings, propagatorSettings_, false, clearNumericalSolution, setIntegratedResult, false,
                        std::chrono::steady_clock::now( ),
                        stateDerivativeModels );

            dynamicsStateDerivative_ = dynamicsSimulator_->getDynamicsStateDerivative( );
            statePostProcessingFunction_ = std::bind(
                        &DynamicsStateDerivativeModel< TimeType, StateScalarType >::postProcessStateAndVariationalEquations,
                        dynamicsStateDerivative_, std::placeholders::_1 );


            // Create variational equations objects.
            variationalEquationsObject_ = std::make_shared< VariationalEquations >(
                        stateDerivativePartials, parametersToEstimate_,
                        dynamicsStateDerivative_->getStateTypeStartIndices( ) );
            dynamicsStateDerivative_->addVariationalEquations( variationalEquationsObject_ );

            // Resize solution of variational equations to 2 (state transition and sensitivity matrices)
            variationalEquationsSolution_.resize( 2 );

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
                            nullptr, nullptr,
                            propagatorSettings_->getConventionalStateSize( ), parameterVectorSize_ );
            }
        }
    }

    //! Destructor
    ~SingleArcVariationalEquationsSolver( ){ }

    //! Function to integrate equations of motion only.
    /*!
     *  Function to integrate equations of motion only (in single arc).  If dynamical
     *  solution is to be processed, the environment is also updated to the new solution.
     *  \param initialStateEstimate Initial state of the equations of motion that is to be used (in same order as in
     *  parametersToEstimate_).
     */
    void integrateDynamicalEquationsOfMotionOnly(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialStateEstimate )
    {
        dynamicsStateDerivative_->setPropagationSettings( std::vector< IntegratedStateType >( ), 1, 0 );
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
    void integrateVariationalAndDynamicalEquations(
            const VectorType& initialStateEstimate, const bool integrateEquationsConcurrently )
    {
        variationalEquationsSolution_[ 0 ].clear( );
        variationalEquationsSolution_[ 1 ].clear( );

        if( integrateEquationsConcurrently )
        {
            // Create initial conditions from new estimate.
            MatrixType initialVariationalState = this->createInitialConditions(
                        dynamicsStateDerivative_->convertFromOutputSolution(
                            initialStateEstimate, integratorSettings_->initialTime_ ) );

            // Integrate variational and state equations.
            dynamicsStateDerivative_->setPropagationSettings( std::vector< IntegratedStateType >( ), 1, 1 );
            dynamicsStateDerivative_->resetFunctionEvaluationCounter( );

            std::map< TimeType, Eigen::VectorXd > dependentVariableHistory;
            std::map< TimeType, MatrixType > rawNumericalSolution;
            std::map< TimeType, double > cumulativeComputationTimeHistory;

            EquationIntegrationInterface< MatrixType, TimeType >::integrateEquations(
                        dynamicsSimulator_->getStateDerivativeFunction( ), rawNumericalSolution,
                        initialVariationalState, integratorSettings_,
                        dynamicsSimulator_->getPropagationTerminationCondition( ),
                        dependentVariableHistory,
                        cumulativeComputationTimeHistory,
                        dynamicsSimulator_->getDependentVariablesFunctions( ),
                        statePostProcessingFunction_,
                        propagatorSettings_->getPrintInterval( ) );

            std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > equationsOfMotionNumericalSolutionRaw;
            std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > equationsOfMotionNumericalSolution;

            utilities::createVectorBlockMatrixHistory(
                        rawNumericalSolution, equationsOfMotionNumericalSolutionRaw,
                        std::make_pair( 0, parameterVectorSize_ ), stateTransitionMatrixSize_ );

            convertNumericalStateSolutionsToOutputSolutions(
                        equationsOfMotionNumericalSolution, equationsOfMotionNumericalSolutionRaw, dynamicsStateDerivative_ );
            dynamicsSimulator_->manuallySetAndProcessRawNumericalEquationsOfMotionSolution(
                        equationsOfMotionNumericalSolution, dependentVariableHistory, dynamicsSimulator_->getSetIntegratedResult( ) );

            // Reset solution for state transition and sensitivity matrices.
            setVariationalEquationsSolution< TimeType, StateScalarType >(
                        rawNumericalSolution, variationalEquationsSolution_,
                        std::make_pair( 0, 0 ), std::make_pair( 0, stateTransitionMatrixSize_ ),
                        stateTransitionMatrixSize_, parameterVectorSize_ );
        }
        else
        {
            dynamicsStateDerivative_->setPropagationSettings( std::vector< IntegratedStateType >( ), 1, 0 );
            dynamicsSimulator_->integrateEquationsOfMotion( initialStateEstimate );

            // Integrate variational equations.
            dynamicsStateDerivative_->setPropagationSettings( { translational_state }, 0, 1 );
            dynamicsStateDerivative_->resetFunctionEvaluationCounter( );

            Eigen::MatrixXd initialVariationalState = this->createInitialVariationalEquationsSolution( );
            std::map< double, Eigen::MatrixXd > rawNumericalSolution;
            std::map< double, Eigen::VectorXd > dependentVariableHistory;
            std::map< double, double > cumulativeComputationTimeHistory;

            EquationIntegrationInterface< Eigen::MatrixXd, double >::integrateEquations(
                        dynamicsSimulator_->getDoubleStateDerivativeFunction( ), rawNumericalSolution, initialVariationalState,
                        variationalOnlyIntegratorSettings_,
                        dynamicsSimulator_->getPropagationTerminationCondition( ),
                        dependentVariableHistory, cumulativeComputationTimeHistory );

            setVariationalEquationsSolution< double, double >(
                        rawNumericalSolution, variationalEquationsSolution_, std::make_pair( 0, 0 ),
                        std::make_pair( 0, stateTransitionMatrixSize_ ),
                        stateTransitionMatrixSize_, parameterVectorSize_ );

        }

        // Reset solution for state transition and sensitivity matrices.
        resetVariationalEquationsInterpolators( );
    }

    //! Function to return the numerical solution history of numerically integrated variational equations.
    /*!
     *  Function to return the numerical solution history of numerically integrated variational equations.
     *  \return Vector of mapa of state transition matrix history (first vector entry)
     *  and sensitivity matrix history (second vector entry)
     */
    std::vector< std::map< double, Eigen::MatrixXd > >& getNumericalVariationalEquationsSolution( )
    {
        return variationalEquationsSolution_;
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
        propagatorSettings_->resetInitialStates(
                    estimatable_parameters::getInitialStateVectorOfBodiesToEstimate( parametersToEstimate_ ) );

        dynamicsStateDerivative_->template updateStateDerivativeModelSettings(
                    propagatorSettings_->getInitialStates( ) );

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
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
                stateTransitionMatrixInterpolator;
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
                sensitivityMatrixInterpolator;
        createStateTransitionAndSensitivityMatrixInterpolator(
                    stateTransitionMatrixInterpolator, sensitivityMatrixInterpolator, variationalEquationsSolution_,
                    this->clearNumericalSolution_ );

        // Create (if non-existent) or reset state transition matrix interface
        if( stateTransitionInterface_ == nullptr )
        {
            stateTransitionInterface_ = std::make_shared< SingleArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                        stateTransitionMatrixInterpolator, sensitivityMatrixInterpolator,
                        propagatorSettings_->getConventionalStateSize( ), parameterVectorSize_ );
        }
        else
        {
            std::dynamic_pointer_cast< SingleArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                        stateTransitionInterface_ )->updateMatrixInterpolators(
                        stateTransitionMatrixInterpolator, sensitivityMatrixInterpolator );
        }
    }

    //! Object used for numerically propagating and managing the solution of the equations of motion.
    std::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > dynamicsSimulator_;

    //!  Object that is used to evaluate the variational equations at the given state and time.
    std::shared_ptr< VariationalEquations > variationalEquationsObject_;

    //! Map of history of numerically integrated variational equations.
    /*!
     *  Map of history of numerically integrated variational equations. Key of map denotes time, values are
     *  state transition matrix Phi (first vector entry) and sensitivity matrix S (second vector entry)
     */
    std::vector< std::map< double, Eigen::MatrixXd > > variationalEquationsSolution_;

    std::function< void( Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& ) > statePostProcessingFunction_;


    //! Settings for numerical integrator of combined propagation of variational equations and equations of motion.
    std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings_;

    //! Settings for propagation of equations of motion.
    std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > propagatorSettings_;

    //! Settings for numerical integrator when integrating only variational equations.
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings_;

    //! Object used to compute the full state derivative in equations of motion and variational equations.
    /*!
     *  Object used to compute the full state derivative in equations of motion and variational equations,
     *  including relevant updates of environment from current state and time. Object may be used for
     *  either full or separate propagation of equations.
     */
    std::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > dynamicsStateDerivative_;
};


//! Function to transfer the initial multi-arc states from propagator settings to associated initial state estimation parameters.
/*!
 *  Function to transfer the initial multi-arc states from propagator settings to associated initial state estimation parameters.
 *  \param parametersToEstimate Full set of estimated parameters to which the initial states are to be transferred
 *  \param propagatorSettings Multi-arc propagator settings from which the initial states are to be taken
 */
template< typename StateScalarType = double >
void setPropagatorSettingsMultiArcStatesInEstimatedDynamicalParameters(
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > >  parametersToEstimate,
        const std::shared_ptr< MultiArcPropagatorSettings< StateScalarType > > propagatorSettings )
{
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > StateType;
    typedef std::map< std::string, std::shared_ptr< estimatable_parameters::EstimatableParameter< StateType > > >
            ArcWiseParameterList;

    // Get list of estimated bodies
    ArcWiseParameterList estimatedBodies = estimatable_parameters::getListOfBodiesWithTranslationalMultiArcStateToEstimate(
                parametersToEstimate );
    std::vector< std::string > bodiesWithPropagatedTranslationalState =
            utilities::createVectorFromMapKeys( estimatedBodies );

    // Iterate over each arc and set initial state.
    std::map< std::string, StateType > arcInitialTranslationalStates;
    for( int arc = 0; arc < propagatorSettings->getNmberOfArcs( ); arc++ )
    {
        // Check type of dynamics
        switch( propagatorSettings->getSingleArcSettings( ).at( arc )->getStateType( ) )
        {
        case translational_state:
        {
            std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > translationalPropagatorSettings =
                    std::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType > >(
                        propagatorSettings->getSingleArcSettings( ).at( arc ) );

            // Iterate over bodies and set initial state
            for( unsigned int i = 0; i < translationalPropagatorSettings->bodiesToIntegrate_.size( ); i++ )
            {
                if( arc == 0 )
                {
                    arcInitialTranslationalStates[ translationalPropagatorSettings->bodiesToIntegrate_.at( i ) ] =
                            StateType( 6 * propagatorSettings->getNmberOfArcs( ) );
                }
                arcInitialTranslationalStates[ translationalPropagatorSettings->bodiesToIntegrate_.at( i ) ].segment( arc * 6, 6 ) =
                        translationalPropagatorSettings->getInitialStates( ).segment( i * 6, 6 );
            }
            break;
        }
        default:
            std::string errorMessage = "Error, cannot yet make parameters and multi-arc propagator settings consistent for " +
                    std::to_string( propagatorSettings->getSingleArcSettings( ).at( arc )->getStateType( ) );
            throw std::runtime_error( errorMessage );
        }
    }

    // Set information in estimation objects
    for( unsigned int i = 0; i < bodiesWithPropagatedTranslationalState.size( ); i++ )
    {
        estimatedBodies.at( bodiesWithPropagatedTranslationalState.at( i ) )->setParameterValue(
                    arcInitialTranslationalStates.at( bodiesWithPropagatedTranslationalState.at( i ) ) );
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
class MultiArcVariationalEquationsSolver: public VariationalEquationsSolver< StateScalarType, TimeType >
{
public:

    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > MatrixType;
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > VectorType;

    using VariationalEquationsSolver< StateScalarType, TimeType >::parametersToEstimate_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::bodyMap_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::stateTransitionMatrixSize_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::parameterVectorSize_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::stateTransitionInterface_;

    //! Constructor
    /*!
     *  Constructor, sets up object for automatic evaluation and numerical integration of variational equations and equations of motion.
     *  \param bodyMap Map of bodies (with names) of all bodies in integration.
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
            const simulation_setup::NamedBodyMap& bodyMap,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
            const std::vector< double > arcStartTimes,
            const bool integrateDynamicalAndVariationalEquationsConcurrently = true,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings =
            std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ),
            const bool clearNumericalSolution = true,
            const bool integrateEquationsOnCreation = false,
            const bool resetMultiArcDynamicsAfterPropagation = true ):
        VariationalEquationsSolver< StateScalarType, TimeType >(
            bodyMap, parametersToEstimate, clearNumericalSolution ),
        propagatorSettings_( std::dynamic_pointer_cast< MultiArcPropagatorSettings< StateScalarType > >( propagatorSettings ) ),
        resetMultiArcDynamicsAfterPropagation_( resetMultiArcDynamicsAfterPropagation )
    {
        if(  std::dynamic_pointer_cast< MultiArcPropagatorSettings< StateScalarType > >( propagatorSettings ) == nullptr )
        {
            throw std::runtime_error( "Error when making multi-arc variational equartions solver, input is single-arc" );
        }
        checkMultiArcPropagatorSettingsAndParameterEstimationConsistency(
                    propagatorSettings_, parametersToEstimate, arcStartTimes );

        parameterVectorSize_ = estimatable_parameters::getSingleArcParameterSetSize( parametersToEstimate );

        stateTransitionMatrixSize_ -= ( parametersToEstimate->getParameterSetSize( ) -
                                        estimatable_parameters::getSingleArcParameterSetSize( parametersToEstimate ) );

        dynamicsSimulator_ =  std::make_shared< MultiArcDynamicsSimulator< StateScalarType, TimeType > >(
                    bodyMap, integratorSettings, propagatorSettings, arcStartTimes,
                    false, clearNumericalSolution, resetMultiArcDynamicsAfterPropagation_ );


        std::vector< std::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > > singleArcDynamicsSimulators =
                dynamicsSimulator_->getSingleArcDynamicsSimulators( );

        if( arcStartTimes.size( ) != singleArcDynamicsSimulators.size( ) )
        {
            throw std::runtime_error( "Error when making multi-arc variational equartions solver, input is inconsistent" );
        }

        for( unsigned int i = 0; i < singleArcDynamicsSimulators.size( ); i++ )
        {
            dynamicsStateDerivatives_.push_back( singleArcDynamicsSimulators.at( i )->getDynamicsStateDerivative( ) );
            // Create variational equations objects.
            std::map< IntegratedStateType, orbit_determination::StateDerivativePartialsMap > stateDerivativePartials =
                    simulation_setup::createStateDerivativePartials< StateScalarType, TimeType >(
                        dynamicsStateDerivatives_.at( i )->getStateDerivativeModels( ), bodyMap, parametersToEstimate );
            std::shared_ptr< VariationalEquations > variationalEquationsObject_ =
                    std::make_shared< VariationalEquations >(
                        stateDerivativePartials, parametersToEstimate_, dynamicsStateDerivatives_.at( i )->getStateTypeStartIndices( ) );

            dynamicsStateDerivatives_.at( i )->addVariationalEquations( variationalEquationsObject_ );
            arcStartTimes_.push_back( arcStartTimes.at( i ) );
        }

        numberOfArcs_ = dynamicsStateDerivatives_.size( );
        // Resize solution of variational equations to 2 (state transition and sensitivity matrices)
        variationalEquationsSolution_.resize( numberOfArcs_ );
        for( int i = 0; i < numberOfArcs_; i++ )
        {
            variationalEquationsSolution_[ i ].resize( 2 );
        }

        // Integrate variational equations from initial state estimate.
        if( integrateEquationsOnCreation )
        {
            if( integrateDynamicalAndVariationalEquationsConcurrently )
            {
                integrateVariationalAndDynamicalEquations( propagatorSettings_->getInitialStateList( ) , 1 );
            }
            else
            {
                integrateVariationalAndDynamicalEquations( propagatorSettings_->getInitialStateList( ), 0 );
            }
        }
    }

    //! Destructor
    /*!
     *  Destructor
     */
    ~MultiArcVariationalEquationsSolver( ){ }

    //! Function to integrate equations of motion only.
    /*!
     *  Function to integrate equations of motion only (for all arcs).  If dynamical
     *  solution is to be processed, the environment is also updated to the new solution.
     *  \param initialStateEstimate Initial state of the equations of motion that is to be used (in same order as in
     *  parametersToEstimate_), concatenated for all arcs.
     */
    void integrateDynamicalEquationsOfMotionOnly(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialStateEstimate )
    {
        for( int i = 0; i < numberOfArcs_; i++ )
        {
            dynamicsStateDerivatives_.at( i )->setPropagationSettings( std::vector< IntegratedStateType >( ), 1, 0 );
        }

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
        for( int i = 0; i < numberOfArcs_; i++ )
        {
            dynamicsStateDerivatives_.at( i )->setPropagationSettings( std::vector< IntegratedStateType >( ), 1, 0 );
        }

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
    void integrateVariationalAndDynamicalEquations(
            const VectorType& concatenatedInitialStates, const bool integrateEquationsConcurrently )
    {
        std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > splitInitialState;

        int currentIndex = 0;
        for( unsigned int i = 0; i < dynamicsSimulator_->getSingleArcDynamicsSimulators( ).size( ); i++ )
        {
            int currentSize = dynamicsSimulator_->getSingleArcDynamicsSimulators( ).at( i )->getPropagatorSettings( )->getConventionalStateSize( );

            splitInitialState.push_back( concatenatedInitialStates.block( currentIndex, 0, currentSize, 1 ) );
            currentIndex += currentSize;
        }

        if( currentIndex != concatenatedInitialStates.rows( ) )
        {
            throw std::runtime_error( "Error when doing multi-arc variational equation integration, "
                                      "input state vector size is incompatible with settings." );
        }

        integrateVariationalAndDynamicalEquations( splitInitialState, integrateEquationsConcurrently );
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
    void integrateVariationalAndDynamicalEquations(
            const std::vector< VectorType >& initialStateEstimate, const bool integrateEquationsConcurrently )
    {
        bool updateInitialStates = false;
        std::vector< VectorType > arcInitialStates;

        // Retrieve single-arc dynamics simulator objects
        std::vector< std::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > > singleArcDynamicsSimulators =
                dynamicsSimulator_->getSingleArcDynamicsSimulators( );

        // Clear solution maps for variational equations
        for( int i = 0; i < numberOfArcs_; i++ )
        {
            variationalEquationsSolution_[ i ][ 0 ].clear( );
            variationalEquationsSolution_[ i ][ 1 ].clear( );
        }

        // Propagate variational equations and equations of motion concurrently
        if( integrateEquationsConcurrently )
        {
            // Allocate maps that stored numerical solution for equations of motion
            std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >
                    currentEquationsOfMotionNumericalSolutionsRaw;
            std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >
                    equationsOfMotionNumericalSolutions;
            std::vector< std::map< TimeType, Eigen::Matrix< double, Eigen::Dynamic, 1 > > >
                    dependentVariableHistorySolutions;
            std::vector< std::map< TimeType, double > > cumulativeComputationTimeHistorySolutions;

            equationsOfMotionNumericalSolutions.resize( numberOfArcs_ );

            dependentVariableHistorySolutions.resize( numberOfArcs_ );
            cumulativeComputationTimeHistorySolutions.resize( numberOfArcs_ );

            // Integrate equations for all arcs.
            for( int i = 0; i < numberOfArcs_; i++ )
            {
                // Retrieve integrator settings, and ensure correct initial time.
                std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings =
                        singleArcDynamicsSimulators.at( i )->getIntegratorSettings( );
                integratorSettings->initialTime_ = arcStartTimes_.at( i );

                // Set state derivative model to propagate both variational equations and equations of motion
                singleArcDynamicsSimulators.at( i )->getDynamicsStateDerivative( )->setPropagationSettings(
                            std::vector< IntegratedStateType >( ), 1, 1 );


                // Get arc initial state. If initial state is NaN, this signals that the initial state is to be taken from
                // previous arc
                VectorType currentArcInitialState;

                if( ( i == 0 ) || ( !linear_algebra::doesMatrixHaveNanEntries( initialStateEstimate.at( i ) ) ) )
                {
                    currentArcInitialState = initialStateEstimate.at( i );
                }
                else
                {
                    currentArcInitialState = getArcInitialStateFromPreviousArcResult(
                                equationsOfMotionNumericalSolutions.at( i - 1 ), arcStartTimes_.at( i ) );
                    updateInitialStates = true;
                }
                arcInitialStates.push_back( currentArcInitialState );

                // Update state derivative model to (possible) update in state.
                singleArcDynamicsSimulators.at( i )->getDynamicsStateDerivative( )->
                        template updateStateDerivativeModelSettings( currentArcInitialState );

                // Create initial state for combined variational/equations of motion.
                MatrixType initialVariationalState = this->createInitialConditions(
                            currentArcInitialState );


                // Integrate variational and state equations.
                dynamicsSimulator_->getDynamicsStateDerivative( ).at( i )->resetFunctionEvaluationCounter( );
                std::map< TimeType, MatrixType > rawNumericalSolution;
                EquationIntegrationInterface< MatrixType, TimeType >::integrateEquations(
                            singleArcDynamicsSimulators.at( i )->getStateDerivativeFunction( ),
                            rawNumericalSolution,
                            initialVariationalState, integratorSettings,
                            singleArcDynamicsSimulators.at( i )->getPropagationTerminationCondition( ),
                            dependentVariableHistorySolutions.at( i ),
                            cumulativeComputationTimeHistorySolutions.at( i ),
                            singleArcDynamicsSimulators.at( i )->getDependentVariablesFunctions( ),
                            std::bind(
                                &DynamicsStateDerivativeModel< TimeType, StateScalarType >::postProcessStateAndVariationalEquations,
                                singleArcDynamicsSimulators.at( i )->getDynamicsStateDerivative( ), std::placeholders::_1 ) );

                // Extract solution of equations of motion.
                utilities::createVectorBlockMatrixHistory(
                            rawNumericalSolution, currentEquationsOfMotionNumericalSolutionsRaw,
                            std::make_pair( 0, parameterVectorSize_ ), stateTransitionMatrixSize_ );


                // Transform equations of motion solution to output formulation
                convertNumericalStateSolutionsToOutputSolutions(
                            equationsOfMotionNumericalSolutions[ i ], currentEquationsOfMotionNumericalSolutionsRaw,
                            dynamicsStateDerivatives_.at( i ) );
                arcStartTimes_[ i ] = equationsOfMotionNumericalSolutions[ i ].begin( )->first;

                // Save state transition and sensitivity matrix solutions for current arc.
                setVariationalEquationsSolution(
                            rawNumericalSolution, variationalEquationsSolution_[ i ],
                            std::make_pair( 0, 0 ), std::make_pair( 0, stateTransitionMatrixSize_ ),
                            stateTransitionMatrixSize_, parameterVectorSize_ );

            }

            // Process numerical solution of equations of motion
            dynamicsSimulator_->manuallySetAndProcessRawNumericalEquationsOfMotionSolution(
                        equationsOfMotionNumericalSolutions, dependentVariableHistorySolutions,
                        resetMultiArcDynamicsAfterPropagation_ );
            equationsOfMotionNumericalSolutions.clear( );

            if( updateInitialStates )
            {
                propagatorSettings_->resetInitialStatesList( arcInitialStates );
                setPropagatorSettingsMultiArcStatesInEstimatedDynamicalParameters(
                            parametersToEstimate_, propagatorSettings_ );
            }
        }
        else
        {
            // Integrate dynamics for each arc
            for( int i = 0; i < numberOfArcs_; i++ )
            {
                // Get arc initial state. If initial state is NaN, this signals that the initial state is to be taken from
                // previous arc
                if( ( i == 0 ) || ( !linear_algebra::doesMatrixHaveNanEntries( initialStateEstimate.at( i ) ) ) )
                {
                    throw std::runtime_error( "Error, arc information transferral not yet supported for separate "
                                              "dynamics and variational euations propagation." );
                    updateInitialStates = true;
                }
                else
                {
                    arcInitialStates.push_back( initialStateEstimate.at( i ) );
                }

                // Update state derivative model to (possible) update in state.
                singleArcDynamicsSimulators.at( i )->getDynamicsStateDerivative( )->
                        template updateStateDerivativeModelSettings( arcInitialStates.at( i ) );
                singleArcDynamicsSimulators.at( i )->getDynamicsStateDerivative( )->setPropagationSettings(
                            std::vector< IntegratedStateType >( ), 1, 0 );
            }

            dynamicsSimulator_->integrateEquationsOfMotion( arcInitialStates );
            arcStartTimes_ = dynamicsSimulator_->getArcStartTimes( );

            std::map< TimeType, MatrixType > rawNumericalSolutions;
            std::map< TimeType, Eigen::Matrix< double, Eigen::Dynamic, 1 > > dummyDependentVariableHistorySolution;
            std::map< TimeType, double > dummyCumulativeComputationTimeHistorySolution;

            // Integrate variational equarions for each arc
            for( int i = 0; i < numberOfArcs_; i++ )
            {
                // Propagate only variational equations
                singleArcDynamicsSimulators.at( i )->getDynamicsStateDerivative( )->setPropagationSettings(
                { translational_state }, 0, 1 );

                // Get initial state for variational equations (single arc)
                MatrixType initialVariationalState = this->createInitialVariationalEquationsSolution( ).
                        template cast< StateScalarType >( );

                // Integrate variational equations for current arc
                dynamicsSimulator_->getDynamicsStateDerivative( ).at( i )->resetFunctionEvaluationCounter( );
                EquationIntegrationInterface< MatrixType, TimeType >::integrateEquations(
                            singleArcDynamicsSimulators.at( i )->getStateDerivativeFunction( ),
                            rawNumericalSolutions, initialVariationalState,
                            singleArcDynamicsSimulators.at( i )->getIntegratorSettings( ),
                            singleArcDynamicsSimulators.at( i )->getPropagationTerminationCondition( ),
                            dummyDependentVariableHistorySolution, dummyCumulativeComputationTimeHistorySolution );

                // Save state transition and sensitivity matrix solutions for current arc.
                setVariationalEquationsSolution(
                            rawNumericalSolutions, variationalEquationsSolution_[ i ],
                            std::make_pair( 0, 0 ), std::make_pair( 0, stateTransitionMatrixSize_ ),
                            stateTransitionMatrixSize_, parameterVectorSize_ );
                rawNumericalSolutions.clear( );
            }

        }

        if( updateInitialStates )
        {
            propagatorSettings_->resetInitialStatesList( arcInitialStates );
            setPropagatorSettingsMultiArcStatesInEstimatedDynamicalParameters(
                        parametersToEstimate_, propagatorSettings_ );
        }

        // Reset solution for state transition and sensitivity matrices.
        resetVariationalEquationsInterpolators( );

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
        propagatorSettings_->resetInitialStates(
                    estimatable_parameters::getInitialStateVectorOfBodiesToEstimate( parametersToEstimate_ ) );


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
        return variationalEquationsSolution_;
    }

    //! Function to return list of start times of each arc. NOTE: This list is updated after every propagation.
    /*!
     * Function to return list of start times of each arc. NOTE: This list is updated after every propagation.
     * \return List of start times of each arc. NOTE: This list is updated after every propagation.
     */
    std::vector< double > getArcStartTimes( )
    {
        return arcStartTimes_;
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
        stateTransitionMatrixInterpolators.resize( variationalEquationsSolution_.size( ) );
        sensitivityMatrixInterpolators.resize( variationalEquationsSolution_.size( ) );

        // Create interpolators.
        for( unsigned int i = 0; i < variationalEquationsSolution_.size( ); i++ )
        {
            createStateTransitionAndSensitivityMatrixInterpolator(
                        stateTransitionMatrixInterpolators[ i ],
                        sensitivityMatrixInterpolators[ i ],
                        variationalEquationsSolution_[ i ],
                        this->clearNumericalSolution_ );
        }

        // Create stare transition matrix interface if needed, reset otherwise.
        if( stateTransitionInterface_ == nullptr )
        {
            stateTransitionInterface_ = std::make_shared< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                        stateTransitionMatrixInterpolators, sensitivityMatrixInterpolators,
                        arcStartTimes_,
                        propagatorSettings_->getSingleArcSettings( ).at( 0 )->getConventionalStateSize( ),
                        parametersToEstimate_->getParameterSetSize( ) );
        }
        else
        {
            std::dynamic_pointer_cast< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                        stateTransitionInterface_ )->updateMatrixInterpolators(
                        stateTransitionMatrixInterpolators, sensitivityMatrixInterpolators,
                        arcStartTimes_ );
        }
    }

    //! Object to propagate the dynamics for all arcs.
    std::shared_ptr< MultiArcDynamicsSimulator< StateScalarType, TimeType > > dynamicsSimulator_;

    //! Numerical solution history of integrated variational equations, per arc.
    /*!
     *  Numerical solution history of integrated variational equations, per arc.
     *  Each vector entry contains the results of a single arc, stored in a vector of maps. Inner vector has size two: first entry
     *  is state transition matrix history, second is sensitivity matrix history, both stored as maps. Key of map denotes time,
     *  values are matrices.
     */
    std::vector< std::vector< std::map< double, Eigen::MatrixXd > > > variationalEquationsSolution_;

    //! List of start times of each arc. NOTE: This list is updated after every propagation.
    std::vector< double > arcStartTimes_;


    //! Settings for propagation of equations of motion.
    std::shared_ptr< MultiArcPropagatorSettings< StateScalarType > > propagatorSettings_;

    //! State derivative models for each arc (retrieved from dynamicsSimulator_).
    std::vector< std::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > > dynamicsStateDerivatives_;

    //! Number of arcs over which propagation is to be performed.
    int numberOfArcs_;

    //! Boolean denoting whether to reset the multi-arc dynamics after propagation.
    const bool resetMultiArcDynamicsAfterPropagation_;

};

//! Class to manage and execute the numerical integration of variational equations of a dynamical system in a combination
//! of single and multiple arcs
/*!
 *  Class to manage and execute the numerical integration of variational equations of a dynamical system, in addition
 *  to the dynamics itself, in a combination  of single and multiple arcs. In this class, the governing equations are set once,
 *  but can be re-integrated for different initial conditions using the same instance of the class.
 */
template< typename StateScalarType = double, typename TimeType = double >
class HybridArcVariationalEquationsSolver: public VariationalEquationsSolver< StateScalarType, TimeType >
{
public:

    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > MatrixType;
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > VectorType;

    using VariationalEquationsSolver< StateScalarType, TimeType >::parametersToEstimate_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::bodyMap_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::stateTransitionMatrixSize_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::parameterVectorSize_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::stateTransitionInterface_;

    //! Constructor
    /*!
     *  Constructor, sets up object for automatic evaluation and numerical integration of variational equations and equations of motion.
     *  \param bodyMap Map of bodies (with names) of all bodies in integration.
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
            const simulation_setup::NamedBodyMap& bodyMap,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
            const std::vector< double > arcStartTimes,
            const bool integrateDynamicalAndVariationalEquationsConcurrently = true,
            const bool clearNumericalSolution = true,
            const bool integrateEquationsOnCreation = false ):
        VariationalEquationsSolver< StateScalarType, TimeType >(
            bodyMap, parametersToEstimate, clearNumericalSolution ),
        integratorSettings_( integratorSettings ),
        arcStartTimes_( arcStartTimes )
    {
        // Cast propagator settings to correct type and check validity
        originalPopagatorSettings_ =
                std::dynamic_pointer_cast< HybridArcPropagatorSettings< StateScalarType > >( propagatorSettings );
        if( originalPopagatorSettings_ == nullptr )
        {
            throw std::runtime_error( "Error when making HybridArcVariationalEquationsSolver, input propagation settings are not hybrid arc" );
        }

        // Get input size of single-arc and input multi-arc
        singleArcDynamicsSize_ = originalPopagatorSettings_->getSingleArcPropagatorSettings( )->getConventionalStateSize( );
        originalMultiArcDynamicsSize_ = originalPopagatorSettings_->getMultiArcPropagatorSettings( )->getConventionalStateSize( );
        originalMultiArcDynamicsSingleArcSize_ = originalPopagatorSettings_->getMultiArcPropagatorSettings( )->getConventionalStateSize( ) /
                arcStartTimes.size( );

        // Create propagator settings with the single arc settings included (at the beginning) in each arc
        std::shared_ptr< MultiArcPropagatorSettings< StateScalarType > > extendedMultiArcSettings =
                getExtendedMultiPropagatorSettings(
                    originalPopagatorSettings_->getSingleArcPropagatorSettings( ),
                    originalPopagatorSettings_->getMultiArcPropagatorSettings( ),
                    arcStartTimes.size( ) );
        multiArcDynamicsSize_ = extendedMultiArcSettings->getConventionalStateSize( );
        multiArcDynamicsSingleArcSize_ = extendedMultiArcSettings->getConventionalStateSize( ) / arcStartTimes_.size( );
        propagatorSettings_ = std::make_shared< HybridArcPropagatorSettings< StateScalarType> >(
                    originalPopagatorSettings_->getSingleArcPropagatorSettings( ), extendedMultiArcSettings );

        // Update estimated parameter vector to extended multi-arc settings
        setExtendedMultiArcParameters( arcStartTimes );

        // Create multi-arc solver with original parameter set
        originalMultiArcSolver_ = std::make_shared< MultiArcVariationalEquationsSolver< StateScalarType, TimeType > >(
                    bodyMap, integratorSettings, originalPopagatorSettings_->getMultiArcPropagatorSettings( ),
                    originalMultiArcParametersToEstimate_, arcStartTimes, integrateDynamicalAndVariationalEquationsConcurrently,
                    std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ),
                    false, false, false );

        // Create variational equations solvers for single- and multi-arc
        integratorSettings->initialTime_ = arcStartTimes.at( 0 );
        singleArcSolver_ = std::make_shared< SingleArcVariationalEquationsSolver< StateScalarType, TimeType > >(
                    bodyMap, integratorSettings, propagatorSettings_->getSingleArcPropagatorSettings( ),
                    singleArcParametersToEstimate_, integrateDynamicalAndVariationalEquationsConcurrently,
                    std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ),
                    false, false );
        multiArcSolver_ = std::make_shared< MultiArcVariationalEquationsSolver< StateScalarType, TimeType > >(
                    bodyMap, integratorSettings, extendedMultiArcSettings,
                    multiArcParametersToEstimate_, arcStartTimes, integrateDynamicalAndVariationalEquationsConcurrently,
                    std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ),
                    false, false, false );

        // Create function to retrieve single-arc initial states for extended multi-arc
        std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > singleArcPropagationSettings =
                std::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType > >(
                    propagatorSettings_->getSingleArcPropagatorSettings( ) );
        if( singleArcPropagationSettings == nullptr )
        {
            throw std::runtime_error( "Error when making HybridArcVariationalEquationsSolver, input single arc is not translational" );
        }
        initialStatesFromSingleArcPropagation_ = std::bind(
                    &getInitialStatesOfBodiesFromFrameManager< TimeType, StateScalarType >,
                    singleArcPropagationSettings->bodiesToIntegrate_,
                    singleArcPropagationSettings->centralBodies_,
                    bodyMap, std::placeholders::_1, createFrameManager( bodyMap ) );

        // Propagate dynamical equations if requested
        if( integrateEquationsOnCreation )
        {
            if( integrateDynamicalAndVariationalEquationsConcurrently )
            {
                integrateVariationalAndDynamicalEquations( propagatorSettings_->getInitialStates( ) , 1 );
            }
            else
            {
                integrateVariationalAndDynamicalEquations( propagatorSettings_->getInitialStates( ), 0 );
            }
        }
    }

    //! Destructor
    ~HybridArcVariationalEquationsSolver( ){ }

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
    void integrateVariationalAndDynamicalEquations(
            const VectorType& initialStateEstimate, const bool integrateEquationsConcurrently )
    {

        // Reset initial time and propagate multi-arc equations
        integratorSettings_->initialTime_ = arcStartTimes_.at( 0 );
        singleArcSolver_->integrateVariationalAndDynamicalEquations(
                    initialStateEstimate.block( 0, 0, singleArcDynamicsSize_, 1 ),
                    integrateEquationsConcurrently );


        // Extract single arc state to update multi-arc initial states
        integratorSettings_->initialTime_ = arcStartTimes_.at( 0 );
        resetMultiArcInitialStates(
                    initialStateEstimate.block( singleArcDynamicsSize_, 0, multiArcDynamicsSize_, 1 ) );

        // Reset initial time and propagate single-arc equations
        integratorSettings_->initialTime_ = arcStartTimes_.at( 0 );
        multiArcSolver_->integrateVariationalAndDynamicalEquations(
                    propagatorSettings_->getMultiArcPropagatorSettings( )->getInitialStates( ),
                    integrateEquationsConcurrently );

        copyExtendedMultiArcInitialStatesToOriginalSettins( );

        // Extract multi-arc solution of dynamics, and remove the single arc bodies from the map.
        std::vector< std::map< TimeType, VectorType > > numericalMultiArcSolution  =
                multiArcSolver_->getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );
        std::vector< std::map< TimeType, Eigen::VectorXd > > dependentVariableHistory  =
                multiArcSolver_->getDynamicsSimulator( )->getDependentVariableHistory( );

        removeSingleArcBodiesFromMultiArcSolultion( numericalMultiArcSolution );

        // Reset original multi-arc bodies' dynamics
        originalMultiArcSolver_->getDynamicsSimulator( )->manuallySetAndProcessRawNumericalEquationsOfMotionSolution(
                    numericalMultiArcSolution, dependentVariableHistory, true );

        // Create state transition matrix if not yet created.
        if( stateTransitionInterface_ == nullptr )
        {
            if( std::dynamic_pointer_cast< SingleArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                        singleArcSolver_->getStateTransitionMatrixInterface( ) ) == nullptr )
            {
                throw std::runtime_error( "Error when making hybrid state transition/sensitivity interface, single-arc input is nullptr" );
            }

            if( std::dynamic_pointer_cast< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                        multiArcSolver_->getStateTransitionMatrixInterface( ) ) == nullptr )
            {
                throw std::runtime_error( "Error when making hybrid state transition/sensitivity interface, multi-arc input is nullptr" );
            }

            stateTransitionInterface_ = std::make_shared< HybridArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                        std::dynamic_pointer_cast< SingleArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                            singleArcSolver_->getStateTransitionMatrixInterface( ) ),
                        std::dynamic_pointer_cast< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface >(
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
    void integrateDynamicalEquationsOfMotionOnly(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialStateEstimate )
    {
        // Reset initial time and propagate multi-arc equations
        integratorSettings_->initialTime_ = arcStartTimes_.at( 0 );
        singleArcSolver_->integrateDynamicalEquationsOfMotionOnly(
                    initialStateEstimate.block( 0, 0, singleArcDynamicsSize_, 1 ) );

        // Extract single arc state to update multi-arc initial states
        integratorSettings_->initialTime_ = arcStartTimes_.at( 0 );
        resetMultiArcInitialStates(
                    initialStateEstimate.block( singleArcDynamicsSize_, 0, multiArcDynamicsSize_, 1 ) );

        // Reset initial time and propagate single-arc equations
        integratorSettings_->initialTime_ = arcStartTimes_.at( 0 );
        multiArcSolver_->integrateDynamicalEquationsOfMotionOnly(
                    propagatorSettings_->getMultiArcPropagatorSettings( )->getInitialStates( ) );

        copyExtendedMultiArcInitialStatesToOriginalSettins( );

        // Extract multi-arc solution of dynamics, and remove the single arc bodies from the map.
        std::vector< std::map< TimeType, VectorType > > numericalMultiArcSolution  =
                multiArcSolver_->getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );
        std::vector< std::map< TimeType, Eigen::VectorXd > > dependentVariableHistory  =
                multiArcSolver_->getDynamicsSimulator( )->getDependentVariableHistory( );
        removeSingleArcBodiesFromMultiArcSolultion( numericalMultiArcSolution );

        // Reset original multi-arc bodies' dynamics
        originalMultiArcSolver_->getDynamicsSimulator( )->manuallySetAndProcessRawNumericalEquationsOfMotionSolution(
                    numericalMultiArcSolution, dependentVariableHistory, true );

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
        originalPopagatorSettings_->resetInitialStates(
                    estimatable_parameters::getInitialStateVectorOfBodiesToEstimate( parametersToEstimate_ ) );

        propagatorSettings_->getSingleArcPropagatorSettings( )->resetInitialStates(
                    newParameterEstimate.segment( 0, singleArcDynamicsSize_ ) );
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > totalMultiArcInitialState =
                propagatorSettings_->getMultiArcPropagatorSettings( )->getInitialStates( );

        for( unsigned int i = 0; i < arcStartTimes_.size( ); i++ )
        {
            totalMultiArcInitialState.segment(
                        i * multiArcDynamicsSingleArcSize_ + singleArcDynamicsSize_,
                        originalMultiArcDynamicsSingleArcSize_ ) =
                    newParameterEstimate.segment(
                        singleArcDynamicsSize_ + i * originalMultiArcDynamicsSingleArcSize_, originalMultiArcDynamicsSingleArcSize_  );

        }
        propagatorSettings_->getMultiArcPropagatorSettings( )->resetInitialStates( totalMultiArcInitialState );
        propagatorSettings_->setInitialStatesFromConstituents( );


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
    std::shared_ptr< HybridArcPropagatorSettings< StateScalarType > > getPropagatorSettings( )
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

        std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
                Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > >
                singleArcParameters = parametersToEstimate_->getEstimatedSingleArcInitialStateParameters( );
        std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
                Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > >
                originalMultiArcParameters = parametersToEstimate_->getEstimatedMultiArcInitialStateParameters( );

        // Get multi-arc parameters associated with estimated single-arc parameters
        std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
                Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > > extendedMultiArcParameters;
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
    }

    //! Function to reset the initial multi-arc states
    /*!
     * Function to reset the initial multi-arc states
     * \param manualMultiArcStates New multi-arc states
     */
    void resetMultiArcInitialStates(
            const VectorType& manualMultiArcStates )
    {
        // Retrieve full multi-arc initial states, with single-arc bodies not (correctly) set
        std::vector< VectorType > arcInitialStates =
                propagatorSettings_->getMultiArcPropagatorSettings( )->getInitialStateList( );

        // Retrieve single-arc states from ephemerides
        int currentArcSize = 0;
        for( unsigned int i = 0; i < arcInitialStates.size( ); i++ )
        {
            currentArcSize = arcInitialStates[ i ].rows( );
            arcInitialStates[ i ].segment( 0, singleArcDynamicsSize_ ) =
                    initialStatesFromSingleArcPropagation_( arcStartTimes_.at( i ) );

            arcInitialStates[ i ].segment( singleArcDynamicsSize_, currentArcSize - singleArcDynamicsSize_ ) =
                    manualMultiArcStates.segment( i * currentArcSize + singleArcDynamicsSize_, currentArcSize - singleArcDynamicsSize_ );
        }

        // Reset initial multi-arc states in propagator settings and estimated parameters
        propagatorSettings_->getMultiArcPropagatorSettings( )->resetInitialStatesList( arcInitialStates );
        propagatorSettings_->setInitialStatesFromConstituents( );
    }

    //! Function that removes the single-arc body data from propagation results before processing data
    /*!
     * Function that removes the single-arc body data from propagation results before processing data
     * \param numericalMultiArcSolution Full numerical solution of single and multi-arc bodies.
     */
    void removeSingleArcBodiesFromMultiArcSolultion(
            std::vector< std::map< TimeType, VectorType > >& numericalMultiArcSolution )
    {
        // Iterate over all arcs
        for( unsigned int i = 0; i < numericalMultiArcSolution.size( ); i++ )
        {
            // Iterate over all times and remove single-arc bodies from solution
            for( typename std::map< TimeType, VectorType >::iterator mapIterator = numericalMultiArcSolution[ i ].begin( );
                 mapIterator != numericalMultiArcSolution[ i ].end( ); mapIterator++ )
            {
                VectorType fullVector = mapIterator->second;numericalMultiArcSolution[ i ][ mapIterator->first ] =
                        fullVector.segment( singleArcDynamicsSize_, originalMultiArcDynamicsSingleArcSize_ );
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
        for( unsigned int i = 0; i < extendedMultiArcInitialStates.size( ); i++ )
        {
            originalMultiArcInitialStates[ i ] = extendedMultiArcInitialStates.at( i ).segment(
                        singleArcDynamicsSize_, originalMultiArcDynamicsSingleArcSize_ );
        }

        originalPopagatorSettings_->getMultiArcPropagatorSettings( )->resetInitialStatesList( originalMultiArcInitialStates );
        originalPopagatorSettings_->setInitialStatesFromConstituents( );
    }

    //! Object to solve multi-arc variational equations (multi-arc bodies only).
    std::shared_ptr< MultiArcVariationalEquationsSolver< StateScalarType, TimeType > > originalMultiArcSolver_;

    //! Object to solve single-arc variational equations.
    std::shared_ptr< SingleArcVariationalEquationsSolver< StateScalarType, TimeType > > singleArcSolver_;

    //! Object to solve multi-arc variational equations (single- and multi-arc bodies).
    std::shared_ptr< MultiArcVariationalEquationsSolver< StateScalarType, TimeType > > multiArcSolver_;

    //! Propagator settings, with single- and multi-arc separately
    std::shared_ptr< HybridArcPropagatorSettings< StateScalarType > > originalPopagatorSettings_;

    //! Propagator settings, with single-arc bodies added to multi-arc list
    std::shared_ptr< HybridArcPropagatorSettings< StateScalarType > > propagatorSettings_;

    //! Settings to be used for integrator
    std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings_;

    //! Size of estimated single-arc dynamical parameters
    int singleArcDynamicsSize_;

    //! Size of single arc of original estimated multi-arc dynamical parameters
    int originalMultiArcDynamicsSingleArcSize_;

    //! Total size of original estimated multi-arc dynamical parameters
    int originalMultiArcDynamicsSize_;

    //! Size of single arc of extended estimated multi-arc dynamical parameters
    int multiArcDynamicsSize_;

    //! Total size of extended estimated multi-arc dynamical parameters
    int multiArcDynamicsSingleArcSize_;

    //! Estimated parameter set with single-arc dynamical parameters only
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > singleArcParametersToEstimate_ ;

    //! Estimated parameter set with original multi-arc dynamical parameters only
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > originalMultiArcParametersToEstimate_;

    //! Estimated parameter set with extended multi-arc dynamical parameters only
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > multiArcParametersToEstimate_ ;

    //! Times at which arcs for multi-arc solution start
    std::vector< double > arcStartTimes_;

    //! Function that retrieves the single-arc bodies' initial states as a function of time
    /*!
     *  Function that retrieves the single-arc bodies' initial states as a function of time, is used to update the multi-arc
     *  initial states after the single-arc propagation
     */
    std::function< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >( const double ) >initialStatesFromSingleArcPropagation_;

};

extern template class VariationalEquationsSolver< double, double >;
extern template class SingleArcVariationalEquationsSolver< double, double >;
extern template class MultiArcVariationalEquationsSolver< double, double >;
extern template class HybridArcVariationalEquationsSolver< double, double >;

#if( BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS )
extern template class VariationalEquationsSolver< long double, double >;
extern template class VariationalEquationsSolver< double, Time >;
extern template class VariationalEquationsSolver< long double, Time >;

extern template class SingleArcVariationalEquationsSolver< long double, double >;
extern template class SingleArcVariationalEquationsSolver< double, Time >;
extern template class SingleArcVariationalEquationsSolver< long double, Time >;

extern template class MultiArcVariationalEquationsSolver< long double, double >;
extern template class MultiArcVariationalEquationsSolver< double, Time >;
extern template class MultiArcVariationalEquationsSolver< long double, Time >;

extern template class HybridArcVariationalEquationsSolver< long double, double >;
extern template class HybridArcVariationalEquationsSolver< double, Time >;
extern template class HybridArcVariationalEquationsSolver< long double, Time >;

#endif

} // namespace propagators

} // namespace tudat




#endif // TUDAT_VARIATIONALEQUATIONSSOLVER_H
