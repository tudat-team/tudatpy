/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_VARIATIONALEQUATIONSSOLVER_BASE_H
#define TUDAT_VARIATIONALEQUATIONSSOLVER_BASE_H

#include "tudat/basics/utilities.h"

#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/math/interpolators/interpolator.h"
#include "tudat/math/basic/linearAlgebra.h"

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/propagators/stateTransitionMatrixInterface.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/astro/ephemerides/tabulatedEphemeris.h"
#include "tudat/simulation/estimation_setup/createStateDerivativePartials.h"
#include "tudat/simulation/estimation_setup/createEstimatableParameters.h"

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
template< typename StateScalarType = double,
          typename TimeType = double,
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
     *  \param bodies Map of bodies (with names) of all bodies in integration.
     *  \param parametersToEstimate Object containing all parameters that are to be estimated and their current
     *  settings and values.
     *  \param clearNumericalSolution Boolean to determine whether to clear the raw numerical solution member variables
     *  (default true) after propagation and resetting of state transition interface.
     */
    VariationalEquationsSolver(
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
            const bool clearNumericalSolution = 1 ):
        parametersToEstimate_( parametersToEstimate ), bodies_( bodies ),
        stateTransitionMatrixSize_( parametersToEstimate_->getInitialDynamicalStateParameterSize( ) ),
        parameterVectorSize_( parametersToEstimate_->getParameterSetSize( ) ), clearNumericalSolution_( clearNumericalSolution )
    {}

    //! Destructor
    virtual ~VariationalEquationsSolver( ) {}

    //! Pure virtual function to integrate variational equations and equations of motion.
    /*!
     *  Pure virtual function to integrate variational equations and equations of motion, to be implemented in derived
     *  class
     *  \param initialStateEstimate Initial state of the equations of motion that is to be used.
     *  \param integrateEquationsConcurrently Variable determining whether the equations of motion are to be
     *  propagated concurrently with variational equations of motion (if true), or before variational equations (if false).
     */
    virtual void integrateVariationalAndDynamicalEquations( const VectorType& initialStateEstimate,
                                                            const bool integrateEquationsConcurrently ) = 0;

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

    virtual std::shared_ptr< SimulationResults< StateScalarType, TimeType > > getVariationalPropagationResults( ) = 0;

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
            throw std::runtime_error( "Error when getting initial condition for variational equations, sizes are incompatible: " +
                                      std::to_string( stateTransitionMatrixSize_ ) + ", " +
                                      std::to_string( initialStateEstimate.rows( ) ) );
        }

        // Initialize initial conditions to zeros.
        MatrixType varSystemInitialState = MatrixType( stateTransitionMatrixSize_, parameterVectorSize_ + 1 ).setZero( );

        // Set initial state transition matrix to identity
        varSystemInitialState.block( 0, 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_ ).setIdentity( );

        // Set initial body states to current estimate of initial body states.
        varSystemInitialState.block( 0, parameterVectorSize_, stateTransitionMatrixSize_, 1 ) = initialStateEstimate;

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
        Eigen::MatrixXd varSystemInitialState = Eigen::MatrixXd::Zero( stateTransitionMatrixSize_, parameterVectorSize_ );

        // Set initial state transition matrix to identity
        varSystemInitialState.block( 0, 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_ ).setIdentity( );

        return varSystemInitialState;
    }

    //! Object containing all parameters that are to be estimated and their current  settings and values.
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate_;

    //! Map of bodies (with names) of all bodies in integration.
    simulation_setup::SystemOfBodies bodies_;

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
        std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >& numericalIntegrationResult,
        std::vector< std::map< double, Eigen::MatrixXd > >& variationalEquationsSolution,
        const std::pair< int, int > stateTransitionStartIndices,
        const std::pair< int, int > sensitivityStartIndices,
        const int stateTransitionMatrixSize,
        const int parameterSetSize )
{
    variationalEquationsSolution.clear( );
    variationalEquationsSolution.resize( 2 );

    for( typename std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >::iterator integrationIterator =
                 numericalIntegrationResult.begin( );
         integrationIterator != numericalIntegrationResult.end( ); )
    {
        // Set result for state transition matrix in each time step.
        variationalEquationsSolution[ 0 ][ integrationIterator->first ] =
                ( integrationIterator->second.block( stateTransitionStartIndices.first,
                                                     stateTransitionStartIndices.second,
                                                     stateTransitionMatrixSize,
                                                     stateTransitionMatrixSize ) )
                        .template cast< double >( );

        // Set result for sensitivity matrix in each time step.
        variationalEquationsSolution[ 1 ][ integrationIterator->first ] =
                ( integrationIterator->second.block( sensitivityStartIndices.first,
                                                     sensitivityStartIndices.second,
                                                     stateTransitionMatrixSize,
                                                     parameterSetSize - stateTransitionMatrixSize ) )
                        .template cast< double >( );
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
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >& stateTransitionMatrixInterpolator,
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >& sensitivityMatrixInterpolator,
        std::map< double, Eigen::MatrixXd >& stateTransitionSolution,
        std::map< double, Eigen::MatrixXd >& sensitivitySolution,
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
        const std::shared_ptr< SingleArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate )
{
    bool isInputConsistent = 1;

    // Check type of dynamics
    switch( propagatorSettings->getStateType( ) )
    {
        case translational_state: {
            std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType, TimeType > > translationalPropagatorSettings =
                    std::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType, TimeType > >( propagatorSettings );

            // Retrieve estimated and propagated translational states, and check equality.
            std::vector< std::string > propagatedBodies = translationalPropagatorSettings->bodiesToIntegrate_;
            std::vector< std::string > estimatedBodies =
                    estimatable_parameters::getListOfBodiesWithTranslationalStateToEstimate( parametersToEstimate );

            if( static_cast< unsigned int >( translationalPropagatorSettings->getPropagatedStateSize( ) ) != propagatedBodies.size( ) * 6 )
            {
                throw std::runtime_error( "Error when propagating variational equations, tranbbslational state vectors not of size 6." );
            }

            if( propagatedBodies.size( ) != estimatedBodies.size( ) )
            {
                std::string errorMessage = "Error, propagated and estimated body vector sizes are inconsistent " +
                        std::to_string( propagatedBodies.size( ) ) + " " + std::to_string( estimatedBodies.size( ) );
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
                                std::string( propagatedBodies.at( i ) ) + " " + std::string( estimatedBodies.at( i ) );
                        throw std::runtime_error( errorMessage );
                        isInputConsistent = 0;
                    }
                }
            }
            break;
        }
        case rotational_state: {
            std::shared_ptr< RotationalStatePropagatorSettings< StateScalarType, TimeType > > rotationalPropagatorSettings =
                    std::dynamic_pointer_cast< RotationalStatePropagatorSettings< StateScalarType, TimeType > >( propagatorSettings );

            // Retrieve estimated and propagated translational states, and check equality.
            std::vector< std::string > propagatedBodies = rotationalPropagatorSettings->bodiesToIntegrate_;
            std::vector< std::string > estimatedBodies =
                    estimatable_parameters::getListOfBodiesWithRotationalStateToEstimate( parametersToEstimate );
            if( propagatedBodies.size( ) != estimatedBodies.size( ) )
            {
                std::string errorMessage = "Error, propagated and estimated body vector sizes are inconsistent " +
                        std::to_string( propagatedBodies.size( ) ) + " " + std::to_string( estimatedBodies.size( ) );
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
                                std::string( propagatedBodies.at( i ) ) + " " + std::string( estimatedBodies.at( i ) );
                        throw std::runtime_error( errorMessage );
                        isInputConsistent = 0;
                    }
                }
            }
            break;
        }
        case body_mass_state: {
            std::shared_ptr< MassPropagatorSettings< StateScalarType, TimeType > > massPropagatorSettings =
                    std::dynamic_pointer_cast< MassPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings );

            // Retrieve estimated and propagated translational states, and check equality.
            std::vector< std::string > propagatedBodies = massPropagatorSettings->bodiesWithMassToPropagate_;
            std::vector< std::string > estimatedBodies =
                    estimatable_parameters::getListOfBodiesWithMassStateToEstimate( parametersToEstimate );
            if( propagatedBodies.size( ) != estimatedBodies.size( ) )
            {
                std::string errorMessage = "Error, propagated and estimated body vector sizes are inconsistent " +
                        std::to_string( propagatedBodies.size( ) ) + " " + std::to_string( estimatedBodies.size( ) );
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
                                std::string( propagatedBodies.at( i ) ) + " " + std::string( estimatedBodies.at( i ) );
                        throw std::runtime_error( errorMessage );
                        isInputConsistent = 0;
                    }
                }
            }
            break;
        }
        case hybrid: {
            std::shared_ptr< MultiTypePropagatorSettings< StateScalarType, TimeType > > multiTypePropagatorSettings =
                    std::dynamic_pointer_cast< MultiTypePropagatorSettings< StateScalarType, TimeType > >( propagatorSettings );
            isInputConsistent = true;

            for( auto settingIterator = multiTypePropagatorSettings->propagatorSettingsMap_.begin( );
                 settingIterator != multiTypePropagatorSettings->propagatorSettingsMap_.end( );
                 settingIterator++ )
            {
                for( unsigned int i = 0; i < settingIterator->second.size( ); i++ )
                {
                    if( !checkPropagatorSettingsAndParameterEstimationConsistency( settingIterator->second.at( i ), parametersToEstimate ) )
                    {
                        isInputConsistent = false;
                    }
                }
            }

            if( estimatable_parameters::getListOfBodiesWithTranslationalStateToEstimate( parametersToEstimate ).size( ) > 0 &&
                multiTypePropagatorSettings->propagatorSettingsMap_.count( translational_state ) == 0 )
            {
                throw std::runtime_error( "Error, estimating but not propagating translational dynamics" );
                isInputConsistent = false;
            }

            if( estimatable_parameters::getListOfBodiesWithRotationalStateToEstimate( parametersToEstimate ).size( ) > 0 &&
                multiTypePropagatorSettings->propagatorSettingsMap_.count( rotational_state ) == 0 )
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
        const std::shared_ptr< MultiArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
        const std::vector< double > arcStartTimes,
        std::map< int, std::vector< std::string > >& estimatedBodiesPerArc,
        std::map< int, std::map< std::string, int > >& arcIndicesPerBody,
        bool& areEstimatedBodiesDifferentPerArc )
{
    bool isInputConsistent = true;

    // Get list of objets and associated bodies to estimate initial arc-wise translational states
    typedef std::map<
            std::string,
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > >
            ArcWiseParameterList;
    ArcWiseParameterList estimatedBodies =
            estimatable_parameters::getListOfBodiesWithTranslationalMultiArcStateToEstimate( parametersToEstimate );

    int numberEstimatedBodies = estimatedBodies.size( );

    // Check that the arc starting times are provided in correct order.
    for( unsigned int i = 0; i < arcStartTimes.size( ) - 1; i++ )
    {
        if( ( arcStartTimes[ i + 1 ] - arcStartTimes[ i ] ) < 0.0 )
        {
            throw std::runtime_error( "Error, arc start times not provided in increasing order." );
        }
    }

    // Initialising vector keeping track of whether each arc is associated with at least one body whose multi-arc state is to be estimated.
    std::vector< bool > detectedEstimatedStatesPerArc;
    for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
    {
        detectedEstimatedStatesPerArc.push_back( false );
    }

    estimatedBodiesPerArc.clear( );
    arcIndicesPerBody.clear( );
    std::vector< int > counterStateIndicesPerBody;
    for( int i = 0; i < numberEstimatedBodies; i++ )
    {
        counterStateIndicesPerBody.push_back( 0 );
    }

    // Iterate over all parameters and check consistency
    unsigned int counterEstimatedBody = 0;
    for( typename ArcWiseParameterList::const_iterator parameterIterator = estimatedBodies.begin( );
         parameterIterator != estimatedBodies.end( );
         parameterIterator++ )
    {
        // Get arc start times of current parameter
        std::vector< double > parameterArcStartTimes =
                std::dynamic_pointer_cast< estimatable_parameters::ArcWiseInitialTranslationalStateParameter< StateScalarType > >(
                        parameterIterator->second )
                        ->getArcStartTimes( );

        // Check that each arc has at least one body whose state is to be estimated.
        for( unsigned int i = 0; i < parameterArcStartTimes.size( ); i++ )
        {
            bool detectedArc = false;
            int indexDetectedArc = 0;
            for( unsigned int j = indexDetectedArc; j < arcStartTimes.size( ); j++ )
            {
                if( std::fabs( arcStartTimes.at( j ) - parameterArcStartTimes.at( i ) ) <
                    std::max( 4.0 * parameterArcStartTimes.at( i ) * std::numeric_limits< double >::epsilon( ), 1.0E-12 ) )
                {
                    detectedArc = true;
                    indexDetectedArc = j;
                    detectedEstimatedStatesPerArc[ j ] = true;

                    estimatedBodiesPerArc[ indexDetectedArc ].push_back( parameterIterator->first );
                    arcIndicesPerBody[ indexDetectedArc ][ parameterIterator->first ] = counterStateIndicesPerBody[ counterEstimatedBody ];
                    counterStateIndicesPerBody[ counterEstimatedBody ] += 1;
                }
            }

            if( !detectedArc )
            {
                isInputConsistent = false;
                throw std::runtime_error( "Error: arc time for " + parameterIterator->first +
                                          " is incompatible with the vector of "
                                          " arc starting times." );
            }
        }

        counterEstimatedBody += 1;
    }

    for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
    {
        if( !detectedEstimatedStatesPerArc[ i ] )
        {
            isInputConsistent = false;
            throw std::runtime_error( "Error, no multi-arc state to be estimated for arc " + std::to_string( i + 1 ) + " out of " +
                                      std::to_string( arcStartTimes.size( ) ) + "." );
        }
    }

    std::map< IntegratedStateType, std::vector< std::string > > propagatedStateTypes;
    std::map< int, std::vector< std::string > > propagatedBodiesPerArc;

    // Iterate over each arc in propagator settings and check consistency
    for( int arc = 0; arc < propagatorSettings->getNmberOfArcs( ); arc++ )
    {
        // Check type of dynamics
        switch( propagatorSettings->getSingleArcSettings( ).at( arc )->getStateType( ) )
        {
            case translational_state: {
                std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType, TimeType > > translationalPropagatorSettings =
                        std::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType, TimeType > >(
                                propagatorSettings->getSingleArcSettings( ).at( arc ) );

                // Retrieve estimated and propagated translational states, and check equality.
                //            std::vector< std::string > propagatedBodies = translationalPropagatorSettings->bodiesToIntegrate_
                propagatedBodiesPerArc[ arc ] = translationalPropagatorSettings->bodiesToIntegrate_;  // propagatedBodies;
                break;
            }
            default:
                std::string errorMessage = "Error, cannot yet check consistency of multi-arc propagator settings for type " +
                        std::to_string( propagatorSettings->getSingleArcSettings( ).at( arc )->getStateType( ) );
                throw std::runtime_error( errorMessage );
        }
    }

    // Check that propagated and estimated bodies are consistent, for each arc.
    for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
    {
        if( estimatedBodiesPerArc.at( i ).size( ) != propagatedBodiesPerArc.at( i ).size( ) )
        {
            isInputConsistent = false;
            throw std::runtime_error( "Error, for arc " + std::to_string( i + 1 ) + " out of " + std::to_string( arcStartTimes.size( ) ) +
                                      ", number of propagated bodies inconsistent with number of estimated bodies." );
        }
        for( unsigned int j = 0; j < estimatedBodiesPerArc.at( i ).size( ); j++ )
        {
            auto currentListToTest = propagatedBodiesPerArc.at( i );
            auto itr = std::find( currentListToTest.begin( ), currentListToTest.end( ), estimatedBodiesPerArc.at( i ).at( j ) );
            if( itr == currentListToTest.end( ) )
            {
                isInputConsistent = false;
                std::string currentPropagatedBodies = "";
                for( unsigned int k = 0; k < currentListToTest.size( ); k++ )
                {
                    currentPropagatedBodies += std::to_string( k ) + ":" + currentListToTest.at( k ) + "; ";
                }
                throw std::runtime_error( "Error, for arc " + std::to_string( i + 1 ) + " out of " +
                                          std::to_string( arcStartTimes.size( ) ) + ", body " + estimatedBodiesPerArc.at( i ).at( j ) +
                                          " is estimated but not propagated.  Propagated bodies are " + currentPropagatedBodies );
            }
        }
    }

    // Check whether the bodies to be estimated are the same for all arcs.
    areEstimatedBodiesDifferentPerArc = false;
    for( unsigned int i = 1; i < arcStartTimes.size( ); i++ )
    {
        // Check if the number of bodies to be estimated is the same for all arcs.
        if( estimatedBodiesPerArc.at( 0 ).size( ) != estimatedBodiesPerArc.at( i ).size( ) )
        {
            areEstimatedBodiesDifferentPerArc = true;
        }
        else  // Check if the names of the estimates are identical for all arcs.
        {
            for( unsigned int j = 0; j < estimatedBodiesPerArc.at( 0 ).size( ); j++ )
            {
                auto itr = std::find(
                        estimatedBodiesPerArc.at( i ).begin( ), estimatedBodiesPerArc.at( i ).end( ), estimatedBodiesPerArc.at( 0 )[ j ] );
                if( itr == estimatedBodiesPerArc.at( i ).end( ) )
                {
                    areEstimatedBodiesDifferentPerArc = true;
                }
            }
        }
    }

    return isInputConsistent;
}

///! Retrieve parameters to be estimated for each arc (arc-wise parameters might differ from one arc to another).
/*!
 * Retrieve parameters to be estimated for each arc (arc-wise parameters might differ from one arc to another).
 * \param parametersToEstimate Pointer for estimated parameters, provided as input of the whole multi-arc variational equations solver.
 * \param arcWiseParametersToEstimate Vector containing the estimated parameters for each arc (returned by reference).
 * \param estimatedBodiesPerArc list of bodies to be estimated, for each arc.
 */
template< typename StateScalarType = double >
void getParametersToEstimatePerArc(
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
        std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > >& arcWiseParametersToEstimate,
        const std::map< int, std::vector< std::string > >& estimatedBodiesPerArc )
{
    // Get list of objets and associated bodies for initial arc-wise translational states to be estimated.
    typedef std::map<
            std::string,
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > >
            ArcWiseParameterList;
    ArcWiseParameterList estimatedBodies =
            estimatable_parameters::getListOfBodiesWithTranslationalMultiArcStateToEstimate( parametersToEstimate );

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > >
            initialStatesParameters = parametersToEstimate->getEstimatedInitialStateParameters( );

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > > doubleParameters =
            parametersToEstimate->getEstimatedDoubleParameters( );

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > > vectorParameters =
            parametersToEstimate->getEstimatedVectorParameters( );

    for( unsigned int i = 0; i < estimatedBodiesPerArc.size( ); i++ )
    {
        std::vector< std::string > arcWiseBodiesToEstimate = estimatedBodiesPerArc.at( i );

        std::vector<
                std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > >
                arcWiseStatesParameters;

        for( unsigned int j = 0; j < initialStatesParameters.size( ); j++ )
        {
            for( unsigned int body = 0; body < arcWiseBodiesToEstimate.size( ); body++ )
            {
                if( arcWiseBodiesToEstimate[ body ] == initialStatesParameters[ j ]->getParameterName( ).second.first )
                {
                    arcWiseStatesParameters.push_back( initialStatesParameters[ j ] );
                }
            }
        }

        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > arcWiseEstimatableParamatersSet =
                std::make_shared< estimatable_parameters::EstimatableParameterSet< StateScalarType > >(
                        doubleParameters, vectorParameters, arcWiseStatesParameters );

        arcWiseParametersToEstimate.push_back( arcWiseEstimatableParamatersSet );
    }
}

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_VARIATIONALEQUATIONSSOLVER_BASE_H
