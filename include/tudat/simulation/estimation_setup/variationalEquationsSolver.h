/*    Copyright (c) 2010-2019, Delft University of Technology
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


#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>

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
        parametersToEstimate_( parametersToEstimate ),
        bodies_( bodies ),
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
            std::to_string( stateTransitionMatrixSize_ ) + ", " + std::to_string( initialStateEstimate.rows( ) ) );
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
    case translational_state:
    {
        std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType, TimeType  > > translationalPropagatorSettings =
                std::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType, TimeType  > >( propagatorSettings );

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
        std::shared_ptr< RotationalStatePropagatorSettings< StateScalarType, TimeType  > > rotationalPropagatorSettings =
                std::dynamic_pointer_cast< RotationalStatePropagatorSettings< StateScalarType, TimeType  > >( propagatorSettings );

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
    case body_mass_state:
    {
        std::shared_ptr< MassPropagatorSettings< StateScalarType, TimeType  > > massPropagatorSettings =
                std::dynamic_pointer_cast< MassPropagatorSettings< StateScalarType, TimeType  > >( propagatorSettings );

        // Retrieve estimated and propagated translational states, and check equality.
        std::vector< std::string > propagatedBodies = massPropagatorSettings->bodiesWithMassToPropagate_;
        std::vector< std::string > estimatedBodies = estimatable_parameters::getListOfBodiesWithMassStateToEstimate(
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
        std::shared_ptr< MultiTypePropagatorSettings< StateScalarType, TimeType  > > multiTypePropagatorSettings =
                std::dynamic_pointer_cast< MultiTypePropagatorSettings< StateScalarType, TimeType  > >( propagatorSettings );
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
        const std::shared_ptr< MultiArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
        const std::vector< double > arcStartTimes,
        std::map< int, std::vector< std::string > >& estimatedBodiesPerArc,
        std::map< int, std::map< std::string, int > >& arcIndicesPerBody,
        bool& areEstimatedBodiesDifferentPerArc )
{
    bool isInputConsistent = true;

    // Get list of objets and associated bodies to estimate initial arc-wise translational states
    typedef std::map< std::string, std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > > ArcWiseParameterList;
    ArcWiseParameterList estimatedBodies = estimatable_parameters::getListOfBodiesWithTranslationalMultiArcStateToEstimate(
                parametersToEstimate );

    int numberEstimatedBodies = estimatedBodies.size( );

    // Check that the arc starting times are provided in correct order.
    for ( unsigned int i = 0 ; i < arcStartTimes.size( ) - 1 ; i++ )
    {
        if ( ( arcStartTimes[ i + 1 ] - arcStartTimes[ i ] ) < 0.0 )
        {
            throw std::runtime_error( "Error, arc start times not provided in increasing order." );
        }
    }

    // Initialising vector keeping track of whether each arc is associated with at least one body whose multi-arc state is to be estimated.
    std::vector< bool > detectedEstimatedStatesPerArc;
    for ( unsigned int i = 0 ; i < arcStartTimes.size( ) ; i++ )
    {
        detectedEstimatedStatesPerArc.push_back( false );
    }

    estimatedBodiesPerArc.clear( );
    arcIndicesPerBody.clear( );
    std::vector< int > counterStateIndicesPerBody;
    for ( int i = 0 ; i < numberEstimatedBodies ; i++ )
    {
        counterStateIndicesPerBody.push_back( 0 );
    }

    // Iterate over all parameters and check consistency
    unsigned int counterEstimatedBody = 0;
    for( typename ArcWiseParameterList::const_iterator parameterIterator = estimatedBodies.begin( ); parameterIterator !=
         estimatedBodies.end( ); parameterIterator++ )
    {
        // Get arc start times of current parameter
        std::vector< double > parameterArcStartTimes =
                std::dynamic_pointer_cast< estimatable_parameters::
                ArcWiseInitialTranslationalStateParameter< StateScalarType > >(
                    parameterIterator->second )->getArcStartTimes( );

        // Check that each arc has at least one body whose state is to be estimated.
        for ( unsigned int i = 0 ; i < parameterArcStartTimes.size( ) ; i++ )
        {
            bool detectedArc = false;
            int indexDetectedArc = 0;
            for ( unsigned int j = indexDetectedArc ; j < arcStartTimes.size( ) ; j++ )
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

            if ( !detectedArc )
            {
                isInputConsistent = false;
                throw std::runtime_error( "Error: arc time for " + parameterIterator->first + " is incompatible with the vector of "
                                                                                              " arc starting times." );
            }
        }

        counterEstimatedBody += 1;
    }

    for ( unsigned int i = 0 ; i < arcStartTimes.size( ) ; i++ )
    {
        if ( !detectedEstimatedStatesPerArc[ i ] )
        {
            isInputConsistent = false;
            throw std::runtime_error( "Error, no multi-arc state to be estimated for arc " + std::to_string( i + 1 ) + " out of "
                                      + std::to_string( arcStartTimes.size( ) ) + "." );
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
        case translational_state:
        {
            std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType, TimeType > > translationalPropagatorSettings =
                    std::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType, TimeType > >(
                        propagatorSettings->getSingleArcSettings( ).at( arc ) );

            // Retrieve estimated and propagated translational states, and check equality.
//            std::vector< std::string > propagatedBodies = translationalPropagatorSettings->bodiesToIntegrate_
            propagatedBodiesPerArc[ arc ] = translationalPropagatorSettings->bodiesToIntegrate_; //propagatedBodies;
            break;
        }
        default:
            std::string errorMessage = "Error, cannot yet check consistency of multi-arc propagator settings for type " +
                    std::to_string( propagatorSettings->getSingleArcSettings( ).at( arc )->getStateType( ) );
            throw std::runtime_error( errorMessage );
        }
    }

    // Check that propagated and estimated bodies are consistent, for each arc.
    for ( unsigned int i = 0 ; i < arcStartTimes.size( ) ; i++ )
    {
        if ( estimatedBodiesPerArc.at( i ).size( ) != propagatedBodiesPerArc.at( i ).size( ) )
        {
            isInputConsistent = false;
            throw std::runtime_error( "Error, for arc " + std::to_string( i+1 ) + " out of " + std::to_string( arcStartTimes.size( ) )
                                      + ", number of propagated bodies inconsistent with number of estimated bodies." );
        }
        for ( unsigned int j = 0 ; j < estimatedBodiesPerArc.at( i ).size( ) ; j++ )
        {
            auto currentListToTest = propagatedBodiesPerArc.at( i );
            auto itr = std::find( currentListToTest.begin( ), currentListToTest.end( ), estimatedBodiesPerArc.at( i ).at( j ) );
            if ( itr == currentListToTest.end( ) )
            {
                isInputConsistent = false;
                std::string currentPropagatedBodies = "";
                for( unsigned int k = 0; k < currentListToTest.size( ); k++ )
                {
                    currentPropagatedBodies += std::to_string( k ) + ":" + currentListToTest.at( k ) + "; ";
                }
                throw std::runtime_error( "Error, for arc " + std::to_string( i + 1 ) + " out of " + std::to_string( arcStartTimes.size( ) )
                                          + ", body " +  estimatedBodiesPerArc.at( i ).at( j ) + " is estimated but not propagated.  Propagated bodies are " + currentPropagatedBodies );
            }
        }
    }

    // Check whether the bodies to be estimated are the same for all arcs.
    areEstimatedBodiesDifferentPerArc = false;
    for ( unsigned int i = 1 ; i < arcStartTimes.size( ) ; i++ )
    {
        // Check if the number of bodies to be estimated is the same for all arcs.
        if (  estimatedBodiesPerArc.at( 0 ).size( ) != estimatedBodiesPerArc.at( i ).size( ) )
        {
            areEstimatedBodiesDifferentPerArc = true;
        }
        else // Check if the names of the estimates are identical for all arcs.
        {
            for ( unsigned int j = 0 ; j < estimatedBodiesPerArc.at( 0 ).size( ) ; j++ )
            {
                auto itr = std::find( estimatedBodiesPerArc.at( i ).begin( ), estimatedBodiesPerArc.at( i ).end( ),
                                      estimatedBodiesPerArc.at( 0 )[  j ] );
                if ( itr == estimatedBodiesPerArc.at( i ).end( ) )
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
    typedef std::map< std::string, std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > > ArcWiseParameterList;
    ArcWiseParameterList estimatedBodies = estimatable_parameters::getListOfBodiesWithTranslationalMultiArcStateToEstimate(
            parametersToEstimate );

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > >
            initialStatesParameters = parametersToEstimate->getEstimatedInitialStateParameters( );

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > > doubleParameters =
            parametersToEstimate->getEstimatedDoubleParameters( );

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > > vectorParameters =
            parametersToEstimate->getEstimatedVectorParameters( );

    for ( unsigned int i = 0 ; i < estimatedBodiesPerArc.size( ) ; i++ )
    {
        std::vector< std::string > arcWiseBodiesToEstimate = estimatedBodiesPerArc.at( i );

        std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > >
                arcWiseStatesParameters;

        for ( unsigned int j = 0 ; j < initialStatesParameters.size( ) ; j++ )
        {
            for ( unsigned int body = 0 ; body < arcWiseBodiesToEstimate.size( ) ; body++ )
            {
                if ( arcWiseBodiesToEstimate[ body ] == initialStatesParameters[ j ]->getParameterName( ).second.first )
                {
                    arcWiseStatesParameters.push_back( initialStatesParameters[ j ] );
                }
            }
        }

        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > arcWiseEstimatableParamatersSet =
                std::make_shared< estimatable_parameters::EstimatableParameterSet< StateScalarType > >
                        ( doubleParameters, vectorParameters, arcWiseStatesParameters );

        arcWiseParametersToEstimate.push_back( arcWiseEstimatableParamatersSet );

    }
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
            bodies, parametersToEstimate, propagatorSettings != nullptr ?
                propagatorSettings->getOutputSettingsWithCheck( )->getClearNumericalSolutions( ) : false ),
        propagatorSettings_( std::dynamic_pointer_cast< SingleArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ) )
    {
        // Check input consistency
        if( std::dynamic_pointer_cast< SingleArcPropagatorSettings< StateScalarType, TimeType >  >( propagatorSettings ) == nullptr )
        {
            throw std::runtime_error( "Error in variational equations solver, input must be single-arc." );
        }
        else if( !checkPropagatorSettingsAndParameterEstimationConsistency< StateScalarType, TimeType >(
                    propagatorSettings_, parametersToEstimate ) )
        {
            throw std::runtime_error(
                        "Error when making single arc variational equations solver, estimated and propagated bodies are inconsistent." );
        }

        // Create state derivative models
        std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > stateDerivativeModels =
                createStateDerivativeModels( propagatorSettings_, bodies, propagatorSettings_->getInitialTime( ) );

        // Create state derivative partials
        std::map< IntegratedStateType, orbit_determination::StateDerivativePartialsMap >
                stateDerivativePartials =
                simulation_setup::createStateDerivativePartials
                < StateScalarType, TimeType >(
                    getStateDerivativeModelMapFromVector( stateDerivativeModels ), bodies, parametersToEstimate );

        // Create object that propagates the dynamics
        dynamicsSimulator_ = std::make_shared< SingleArcDynamicsSimulator< StateScalarType, TimeType > >(
                    bodies, propagatorSettings_, false,
                    PredefinedSingleArcStateDerivativeModels< StateScalarType, TimeType >( stateDerivativeModels, stateDerivativePartials ) );

        // Create variational equations evaluation objects.
        variationalEquationsObject_ = std::make_shared< VariationalEquations >(
                    stateDerivativePartials, parametersToEstimate_,
                    dynamicsSimulator_->getDynamicsStateDerivative( )->getStateTypeStartIndices( ) );
        dynamicsSimulator_->getDynamicsStateDerivative( )->addVariationalEquations( variationalEquationsObject_ );

        // Create object that will contain and process the propagation results
        variationalPropagationResults_ = std::make_shared< SingleArcVariationalSimulationResults< StateScalarType, TimeType>>(
                dynamicsSimulator_->getSingleArcPropagationResults( ),
                this->stateTransitionMatrixSize_, this->parameterVectorSize_ - this->stateTransitionMatrixSize_ );

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
                        propagatorSettings_->getConventionalStateSize( ), parameterVectorSize_,
                        variationalEquationsObject_->getStatePartialAdditionIndices( ) );
        }
    }

    SingleArcVariationalEquationsSolver(
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
            const bool integrateDynamicalAndVariationalEquationsConcurrently = true,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings
            = std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ),
            const bool clearNumericalSolution = true,
            const bool integrateEquationsOnCreation = true,
            const bool setIntegratedResult = true,
            const bool printDependentVariableData = true,
            const bool setDependentVariablesInterface = false ):
        SingleArcVariationalEquationsSolver( bodies,validateDeprecatedSingleArcSettings(
                                                 integratorSettings, propagatorSettings,
                                                 clearNumericalSolution, setIntegratedResult, false,
                                                 printDependentVariableData, false, setDependentVariablesInterface ), parametersToEstimate,
                                             integrateDynamicalAndVariationalEquationsConcurrently, integrateEquationsOnCreation ){ }

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
        if( integrateEquationsConcurrently )
        {
            // Create initial conditions from new estimate.
            MatrixType initialVariationalState = this->createInitialConditions(
                        dynamicsSimulator_->getDynamicsStateDerivative( )->convertFromOutputSolution(
                            initialStateEstimate, propagatorSettings_->getInitialTime( ) ) );

            // Propagate dynamics and variational equations
            dynamicsSimulator_->integrateEquationsOfMotion( initialVariationalState, variationalPropagationResults_ );
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
    std::vector< std::map< double, Eigen::MatrixXd > > getNumericalVariationalEquationsSolution( )
    {
        std::cerr<<"Warning, use of deprecated single-arc getNumericalVariationalEquationsSolution is not recommended"<<std::endl;
        return std::vector< std::map< double, Eigen::MatrixXd > >( { getStateTransitionMatrixSolution( ), getSensitivityMatrixSolution( ) } );
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
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
                stateTransitionMatrixInterpolator;
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
                sensitivityMatrixInterpolator;

        try
        {
            createStateTransitionAndSensitivityMatrixInterpolator(
                        stateTransitionMatrixInterpolator, sensitivityMatrixInterpolator,
                        variationalPropagationResults_->getStateTransitionSolution( ),
                        variationalPropagationResults_->getSensitivitySolution( ),
                        this->clearNumericalSolution_ );

        }
        catch( const std::exception& caughtException )
        {
            std::cerr << "Error occured when post-processing single-arc variational equation integration results, and creating interpolators, caught error is: " << std::endl << std::endl;
            std::cerr << caughtException.what( ) << std::endl << std::endl;
            std::cerr << "The problem may be that there is an insufficient number of data points (epochs) at which propagation results are produced for one or more arcs. Integrated results are given at" +
                         std::to_string( variationalPropagationResults_->getStateTransitionSolution( ).size( ) ) + " epochs"<< std::endl;
        }


        // Create (if non-existent) or reset state transition matrix interface
        if( stateTransitionInterface_ == nullptr )
        {
            stateTransitionInterface_ = std::make_shared< SingleArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                        stateTransitionMatrixInterpolator, sensitivityMatrixInterpolator,
                        propagatorSettings_->getConventionalStateSize( ), parameterVectorSize_,
                        variationalEquationsObject_->getStatePartialAdditionIndices( ) );
        }
        else
        {
            std::dynamic_pointer_cast< SingleArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                        stateTransitionInterface_ )->updateMatrixInterpolators(
                        stateTransitionMatrixInterpolator, sensitivityMatrixInterpolator,
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


//! Function to transfer the initial multi-arc states from propagator settings to associated initial state estimation parameters.
/*!
 *  Function to transfer the initial multi-arc states from propagator settings to associated initial state estimation parameters.
 *  \param parametersToEstimate Full set of estimated parameters to which the initial states are to be transferred
 *  \param propagatorSettings Multi-arc propagator settings from which the initial states are to be taken
 */
template< typename StateScalarType, typename TimeType >
void setPropagatorSettingsMultiArcStatesInEstimatedDynamicalParameters(
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > >  parametersToEstimate,
        const std::shared_ptr< MultiArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings )
{
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > StateType;
    typedef std::map< std::string, std::shared_ptr< estimatable_parameters::EstimatableParameter< StateType > > > ArcWiseParameterList;

    // Get list of estimated bodies
    ArcWiseParameterList estimatedBodies = estimatable_parameters::getListOfBodiesWithTranslationalMultiArcStateToEstimate(
                parametersToEstimate );
    std::vector< std::string > bodiesWithPropagatedTranslationalState =
            utilities::createVectorFromMapKeys( estimatedBodies );

    std::map< std::string, unsigned int > counterArcPerBody;
    std::map< std::string, StateType > arcInitialTranslationalStates;

    // Iterate over each arc and set initial state.
    std::map< std::string, std::vector< StateType > > arcInitialTranslationalStatesVector;
    for( int arc = 0; arc < propagatorSettings->getNmberOfArcs( ); arc++ )
    {
        // Check type of dynamics
        switch( propagatorSettings->getSingleArcSettings( ).at( arc )->getStateType( ) )
        {
            case translational_state:
            {
                std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType, TimeType > > translationalPropagatorSettings =
                        std::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType, TimeType > >(
                            propagatorSettings->getSingleArcSettings( ).at( arc ) );

                std::vector< std::string > bodiesToIntegrate = translationalPropagatorSettings->bodiesToIntegrate_;

                // Iterate over bodies and set initial state
                for( unsigned int i = 0; i < bodiesToIntegrate.size( ); i++ )
                {
                    if ( counterArcPerBody.count( bodiesToIntegrate.at( i ) ) == 0 )
                    {
                        counterArcPerBody[ bodiesToIntegrate.at( i ) ] = 0;
                        arcInitialTranslationalStatesVector[ bodiesToIntegrate.at( i ) ] = {  };
                    }
                    else
                    {
                        counterArcPerBody.at( bodiesToIntegrate.at( i ) ) += 1;
                    }
                    arcInitialTranslationalStatesVector.at( bodiesToIntegrate.at( i ) ).push_back( translationalPropagatorSettings->getInitialStates( ).segment( i * 6, 6 ) );
                }
                break;
            }
            default:
                std::string errorMessage = "Error, cannot yet make parameters and multi-arc propagator settings consistent for " +
                        std::to_string( propagatorSettings->getSingleArcSettings( ).at( arc )->getStateType( ) );
                throw std::runtime_error( errorMessage );
        }
    }

    for ( auto itr : arcInitialTranslationalStatesVector )
    {
        arcInitialTranslationalStates[ itr.first ] = StateType( 6 * itr.second.size( ) );
        for ( unsigned int k = 0 ; k < itr.second.size( ) ; k++ )
        {
            arcInitialTranslationalStates.at( itr.first ).segment( k*6, 6 ) = itr.second.at( k );
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
            bodies, parametersToEstimate, propagatorSettings != nullptr ?
                propagatorSettings->getOutputSettingsWithCheck( )->getClearNumericalSolutions( ) : false  ),
        propagatorSettings_( propagatorSettings ),
        resetMultiArcDynamicsAfterPropagation_( propagatorSettings != nullptr ?
                propagatorSettings->getOutputSettingsWithCheck( )->getSetIntegratedResult( ) : false )
    {
        if(  std::dynamic_pointer_cast< MultiArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ) == nullptr )
        {
            throw std::runtime_error( "Error when making multi-arc variational equations solver, input is single-arc" );
        }
        checkMultiArcPropagatorSettingsAndParameterEstimationConsistency(
                    propagatorSettings_, parametersToEstimate, propagatorSettings->getArcStartTimes( ), estimatedBodiesPerArc_, arcIndicesPerBody_,
                    areEstimatedBodiesDifferentPerArc_ );

        if ( areEstimatedBodiesDifferentPerArc_ && resetMultiArcDynamicsAfterPropagation_ )
        {
//            throw std::runtime_error( "Error in multi-arc variational equations solver, boolean resetMultiArcDynamicsAfterPropagation should be "
//                                      "set to false when the bodies to be estimated differ from one arc to another." );
        }

        arcWiseParametersToEstimate_.clear( );
        estimatable_parameters::getParametersToEstimatePerArcTest(
                    parametersToEstimate, arcWiseParametersToEstimate_,
                                           propagatorSettings->getArcStartTimes( ),
                                           estimatedBodiesPerArc_, arcIndicesPerBody_ );

        parameterVectorSize_ = 0;
        stateTransitionMatrixSize_  = 0;

        for ( unsigned int arc = 0 ; arc < estimatedBodiesPerArc_.size( ) ; arc++ )
        {
            arcWiseStateTransitionMatrixSize_.push_back( estimatable_parameters::getSingleArcInitialDynamicalStateParameterSetSize( parametersToEstimate, arc ) );
            arcWiseParameterVectorSize_.push_back( estimatable_parameters::getSingleArcParameterSetSize( parametersToEstimate, arc ) );
        }

        dynamicsSimulator_ =  std::make_shared< MultiArcDynamicsSimulator< StateScalarType, TimeType > >(
                    bodies, propagatorSettings, false );

        std::vector< std::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > > singleArcDynamicsSimulators =
                dynamicsSimulator_->getSingleArcDynamicsSimulators( );

        std::vector< std::shared_ptr< SingleArcVariationalSimulationResults< StateScalarType, TimeType > > > singleArcVariationalResults;
        for( unsigned int i = 0; i < dynamicsSimulator_->getSingleArcDynamicsSimulators( ).size( ); i++ )
        {
            singleArcVariationalResults.push_back( std::make_shared< SingleArcVariationalSimulationResults< StateScalarType, TimeType > >(
                    dynamicsSimulator_->getSingleArcDynamicsSimulators( ).at( i )->getSingleArcPropagationResults( ),
                    arcWiseStateTransitionMatrixSize_.at( i ), arcWiseParameterVectorSize_.at( i ) - arcWiseStateTransitionMatrixSize_.at( i ) ) );
        }
        variationalPropagationResults_ = std::make_shared< MultiArcSimulationResults< SingleArcVariationalSimulationResults, StateScalarType, TimeType > >(
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
                    std::make_shared< VariationalEquations >(
                        stateDerivativePartials, arcWiseParametersToEstimate_[ i ], dynamicsStateDerivatives_.at( i )->getStateTypeStartIndices( ), i,
                        arcIndicesPerBody_[ i ] );

            dynamicsStateDerivatives_.at( i )->addVariationalEquations( variationalEquationsObject_ );
        }

        numberOfArcs_ = dynamicsStateDerivatives_.size( );
        // Integrate variational equations from initial state estimate.
        if( integrateEquationsOnCreation )
        {
            integrateVariationalAndDynamicalEquations( propagatorSettings_->getInitialStateList( ) , 1 );
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
        MultiArcVariationalEquationsSolver( bodies, validateDeprecatedMultiArcSettings(
                                        integratorSettings, propagatorSettings, propagationStartTimes,
                                        clearNumericalSolution, resetMultiArcDynamicsAfterPropagation, setDependentVariablesInterface ),
                                            parametersToEstimate,
                                    integrateEquationsOnCreation ){ }


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


    std::shared_ptr< MultiArcInitialStateProvider< StateScalarType > > getInitialStateProvider(
            const std::vector< VectorType >& initialStateEstimate )
    {
        std::vector< std::pair< int, int > > variationalEquationsSize = utilities::mergeVectorsIntoVectorOfPairs(
                arcWiseStateTransitionMatrixSize_, arcWiseParameterVectorSize_ );
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
    void integrateVariationalAndDynamicalEquations(
            const std::vector< VectorType >& initialStateEstimate, const bool integrateEquationsConcurrently )
    {

        // Propagate variational equations and equations of motion concurrently
        if( integrateEquationsConcurrently )
        {
            // Propagate dynamics and variational equations and store results in variationalPropagationResults_ object
            dynamicsSimulator_->template integrateEquationsOfMotion<MultiArcVariationalResults>(
                    variationalPropagationResults_, getInitialStateProvider( initialStateEstimate ));

            // Ensure consistency between parameters and propagator settings
            setPropagatorSettingsMultiArcStatesInEstimatedDynamicalParameters<StateScalarType, TimeType>(
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
        simulation_setup::setInitialStateVectorFromParameterSet< StateScalarType, TimeType >( parametersToEstimate_, propagatorSettings_ );

        for ( int i = 0 ; i < numberOfArcs_ ; i++ )
        {
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >  newParametersValues = propagatorSettings_->getSingleArcSettings( ).at( i )->getInitialStates( );
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > arcWiseParametersValues = arcWiseParametersToEstimate_.at( i )->template getFullParameterValues< StateScalarType >( );
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
        std::cerr<<"Warning, use of deprecated multi-arc getNumericalVariationalEquationsSolution is not recommended"<<std::endl;
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

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > > getArcWiseParametersToEstimate( ) const
    {
        return arcWiseParametersToEstimate_;
    }


   std::shared_ptr< MultiArcVariationalResults > getMultiArcVariationalPropagationResults()
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
                std::cerr << "Error occured when post-processing multi-arc variational equation integration results, and creating interpolators in arc" + std::to_string( i ) + ", caught error is: " << std::endl << std::endl;
                std::cerr << caughtException.what( ) << std::endl << std::endl;
                std::cerr << "The problem may be that there is an insufficient number of data points (epochs) at which propagation results are produced for one or more arcs. Integrated results are given at" +
                             std::to_string( variationalPropagationResults_->getSingleArcResults( ).at( 0 )->getStateTransitionSolution( ).size( ) ) + " epochs"<< std::endl;
            }


        }

        // Create stare transition matrix interface if needed, reset otherwise.
        if( stateTransitionInterface_ == nullptr )
        {
            stateTransitionInterface_ = std::make_shared< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface< StateScalarType > >(
                        stateTransitionMatrixInterpolators, sensitivityMatrixInterpolators,
                        propagatorSettings_->getArcStartTimes( ),
                        arcStartTimesToUse,
                        arcEndTimesToUse,
                        parametersToEstimate_,
                        propagatorSettings_->getSingleArcSettings( ).at( 0 )->getConventionalStateSize( ),
                        parametersToEstimate_->getParameterSetSize( ), getArcWiseStatePartialAdditionIndices( ) );
        }
        else
        {
            std::dynamic_pointer_cast< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface< StateScalarType > >(
                        stateTransitionInterface_ )->updateMatrixInterpolators(
                        stateTransitionMatrixInterpolators, sensitivityMatrixInterpolators,
                        arcStartTimesToUse, arcEndTimesToUse, getArcWiseStatePartialAdditionIndices( ) );
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
            bodies, parametersToEstimate, propagatorSettings != nullptr ?
                propagatorSettings->getOutputSettingsWithCheck( )->getClearNumericalSolutions( ) : false )
    {
        initializeHybridArcVariationalEquationsSolver(
                    bodies, propagatorSettings, integrateEquationsOnCreation );
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
            bodies,  validateDeprecatedHybridArcSettings< StateScalarType, TimeType >(
                integratorSettings,  propagatorSettings,  arcStartTimes, clearNumericalSolution, true, setDependentVariablesInterface ),
            parametersToEstimate, integrateEquationsOnCreation ){ }

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
            bodies, validateDeprecatedHybridArcSettings< StateScalarType, TimeType >(
                singleArcIntegratorSettings,  multiArcIntegratorSettings, propagatorSettings,  arcStartTimes, clearNumericalSolution, true, setDependentVariablesInterface ),
            parametersToEstimate, integrateEquationsOnCreation ){ }

    void initializeHybridArcVariationalEquationsSolver(
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const bool integrateEquationsOnCreation )
    {
        // Cast propagator settings to correct type and check validity
        originalPopagatorSettings_ =
                std::dynamic_pointer_cast< HybridArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings );
        if( originalPopagatorSettings_ == nullptr )
        {
            throw std::runtime_error( "Error when making HybridArcVariationalEquationsSolver, input propagation settings are not hybrid arc" );
        }

        // Retrive arc properties
        singleArcInitialTime_ = originalPopagatorSettings_->getSingleArcPropagatorSettings( )->getInitialTime( );
        int numberOfArcs = originalPopagatorSettings_->getMultiArcPropagatorSettings( )->getNmberOfArcs( );
        arcStartTimes_ = estimatable_parameters::getMultiArcStateEstimationArcStartTimes(
                            parametersToEstimate_, false );

        // Get input size of single-arc and input multi-arc
        singleArcDynamicsSize_ = originalPopagatorSettings_->getSingleArcPropagatorSettings( )->getConventionalStateSize( );
        originalMultiArcDynamicsSize_ = originalPopagatorSettings_->getMultiArcPropagatorSettings( )->getConventionalStateSize( );
        for ( unsigned int i = 0 ; i < arcStartTimes_.size( ) ; i++ )
        {
            originalMultiArcDynamicsSingleArcSize_.push_back(
                    originalPopagatorSettings_->getMultiArcPropagatorSettings( )->getSingleArcSettings( ).at( i )->getConventionalStateSize( ) );
        }

        // Create propagator settings with the single arc settings included (at the beginning) in each arc
        std::shared_ptr< MultiArcPropagatorSettings< StateScalarType, TimeType > > extendedMultiArcSettings =
                getExtendedMultiPropagatorSettings(
                    originalPopagatorSettings_->getSingleArcPropagatorSettings( ),
                    originalPopagatorSettings_->getMultiArcPropagatorSettings( ),
                    numberOfArcs );

        multiArcDynamicsSize_ = extendedMultiArcSettings->getConventionalStateSize( );
        for ( unsigned int i = 0 ; i < arcStartTimes_.size( ) ; i++ )
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
                    bodies, originalPopagatorSettings_->getMultiArcPropagatorSettings( ),
                    originalMultiArcParametersToEstimate_, false );

        // Create variational equations solvers for single- and multi-arc
        singleArcSolver_ = std::make_shared< SingleArcVariationalEquationsSolver< StateScalarType, TimeType > >(
                    bodies, originalPopagatorSettings_->getSingleArcPropagatorSettings( ),
                    singleArcParametersToEstimate_, true, false );

        multiArcSolver_ = std::make_shared< MultiArcVariationalEquationsSolver< StateScalarType, TimeType > >(
                    bodies, extendedMultiArcSettings,
                    multiArcParametersToEstimate_, false );

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

        variationalPropagationResults_ = std::make_shared< HybridArcResults >(
                singleArcSolver_->getSingleArcVariationalPropagationResults( ),
                originalMultiArcSolver_->getMultiArcVariationalPropagationResults( ) );
        initialStatesFromSingleArcPropagation_ = std::bind(
                    &getInitialStatesOfBodiesFromFrameManager< TimeType, StateScalarType >,
                    singleArcPropagationSettings->bodiesToIntegrate_,
                    singleArcPropagationSettings->centralBodies_,
                    bodies, std::placeholders::_1, createFrameManager( bodies.getMap( ) ) );


        // Propagate dynamical equations if requested
        if( integrateEquationsOnCreation )
        {
            integrateVariationalAndDynamicalEquations( propagatorSettings_->getInitialStates( ) , 1 );
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
        // TODO: do process depdendent variables in original multi-arc solver, do not process dependent variables in
        // extended solver. Also add dependent variables to original multi-arc solver.

        // Reset initial time and propagate multi-arc equations
        singleArcSolver_->integrateVariationalAndDynamicalEquations(
                    initialStateEstimate.block( 0, 0, singleArcDynamicsSize_, 1 ),
                    integrateEquationsConcurrently );

        // Extract single arc state to update multi-arc initial states
        resetMultiArcInitialStates(
                    initialStateEstimate.block( singleArcDynamicsSize_, 0, multiArcDynamicsSize_, 1 ) );

        // Reset initial time and propagate single-arc equations
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

        originalMultiArcSolver_->getDynamicsSimulator( )->getMultiArcPropagationResults( )->restartPropagation();

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

            stateTransitionInterface_ = std::make_shared< HybridArcCombinedStateTransitionAndSensitivityMatrixInterface< StateScalarType > >(
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
    void integrateDynamicalEquationsOfMotionOnly(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialStateEstimate )
    {
        // Reset initial time and propagate multi-arc equations
        singleArcSolver_->integrateDynamicalEquationsOfMotionOnly(
                    initialStateEstimate.block( 0, 0, singleArcDynamicsSize_, 1 ) );

        // Extract single arc state to update multi-arc initial states
        resetMultiArcInitialStates(
                    initialStateEstimate.block( singleArcDynamicsSize_, 0, multiArcDynamicsSize_, 1 ) );

        // Reset initial time and propagate single-arc equations
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
        originalMultiArcSolver_->getDynamicsSimulator( )->getMultiArcPropagationResults( )->manuallySetPropagationResults( numericalMultiArcSolution );
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

        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > newInitialStatesOriginalPropagatorSettings = originalPopagatorSettings_->getInitialStates( );

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
                        singleArcDynamicsSize_ + counterOriginalArcWiseIndex /*i * originalMultiArcDynamicsSingleArcSize_.at( i )*/, originalMultiArcDynamicsSingleArcSize_.at( i )  );

            counterFullArcWiseIndex += multiArcDynamicsSingleArcSize_.at( i );
            counterOriginalArcWiseIndex += originalMultiArcDynamicsSingleArcSize_.at( i );

        }
        propagatorSettings_->getMultiArcPropagatorSettings( )->resetInitialStates( totalMultiArcInitialState );
        propagatorSettings_->setInitialStatesFromConstituents( );

        // Reset parameters for arc-wise parameters in both originalMultiArcSolver_ and multiArcSolver_
        for ( unsigned int i = 0 ; i < arcStartTimes_.size( ) ; i++ )
        {
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > newParametersValues = originalPopagatorSettings_->getMultiArcPropagatorSettings( )->getSingleArcSettings( ).at( i )->getInitialStates( );
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > newArcWiseParametersValues = originalMultiArcSolver_->getArcWiseParametersToEstimate( ).at( i )->template getFullParameterValues< StateScalarType >( );
            newArcWiseParametersValues.segment( 0, newParametersValues.size( ) ) = newParametersValues;
            originalMultiArcSolver_->getArcWiseParametersToEstimate( ).at( i )->resetParameterValues( newArcWiseParametersValues );

            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > newFullParametersValues = propagatorSettings_->getMultiArcPropagatorSettings( )->getSingleArcSettings( ).at( i )->getInitialStates( );
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > newFullArcWiseParametersValues = multiArcSolver_->getArcWiseParametersToEstimate( ).at( i )->template getFullParameterValues< StateScalarType >( );
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

    std::shared_ptr< HybridArcResults > getHybridArcVariationalPropagationResults()
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
//        std::cout << "TEST: " << "\n\n";
        estimatable_parameters::printEstimatableParameterEntries( multiArcParametersToEstimate_ );
    }

    //! Function to reset the initial multi-arc states
    /*!
     * Function to reset the initial multi-arc states
     * \param manualMultiArcStates New multi-arc states
     */
    void resetMultiArcInitialStates(
            const VectorType& manualMultiArcStates )
    {
//        std::cout << "manualMultiArcStates: " << manualMultiArcStates.transpose( ) << "\n\n";
//        std::cout << "size manualMultiArcStates: " << manualMultiArcStates.size( ) << "\n\n";

        // Retrieve full multi-arc initial states, with single-arc bodies not (correctly) set
        std::vector< VectorType > arcInitialStates =
                propagatorSettings_->getMultiArcPropagatorSettings( )->getInitialStateList( );

        // Retrieve single-arc states from ephemerides
        int currentArcSize = 0;
        int indexMultiArcState = 0;
        for( unsigned int i = 0; i < arcInitialStates.size( ); i++ )
        {
            currentArcSize = arcInitialStates[ i ].rows( );
//            std::cout << "currentArcSize: " << currentArcSize << "\n\n";
//            std::cout << "singleArcDynamicsSize_" << singleArcDynamicsSize_ << "\n\n";
            arcInitialStates[ i ].segment( 0, singleArcDynamicsSize_ ) =
                    initialStatesFromSingleArcPropagation_( arcStartTimes_.at( i ) );
//            std::cout << "new single arc initial state:" << initialStatesFromSingleArcPropagation_( arcStartTimes_.at( i ) ).transpose( ) << "\n\n";

//            std::cout << "currentArcSize - singleArcDynamicsSize_: " << currentArcSize - singleArcDynamicsSize_ << "\n\n";
//            std::cout << "i * currentArcSize + singleArcDynamicsSize_: " << i * currentArcSize + singleArcDynamicsSize_ << "\n\n";
            arcInitialStates[ i ].segment( singleArcDynamicsSize_, currentArcSize - singleArcDynamicsSize_ ) =
                    manualMultiArcStates.segment( indexMultiArcState /*i * currentArcSize*/ + singleArcDynamicsSize_, currentArcSize - singleArcDynamicsSize_ );
            indexMultiArcState += currentArcSize;
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
//        std::cout << "size numerical multi-arc solution: " << numericalMultiArcSolution.size( ) << "\n\n";
        for( unsigned int i = 0; i < numericalMultiArcSolution.size( ); i++ )
        {
            // Iterate over all times and remove single-arc bodies from solution
            for( typename std::map< TimeType, VectorType >::iterator mapIterator = numericalMultiArcSolution[ i ].begin( );
                 mapIterator != numericalMultiArcSolution[ i ].end( ); mapIterator++ )
            {
                VectorType fullVector = mapIterator->second;numericalMultiArcSolution[ i ][ mapIterator->first ] =
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
            originalMultiArcInitialStates[ i ] = extendedMultiArcInitialStates.at( i ).segment(
                        singleArcDynamicsSize_, originalMultiArcDynamicsSingleArcSize_.at( i ) );
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
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > singleArcParametersToEstimate_ ;

    //! Estimated parameter set with original multi-arc dynamical parameters only
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > originalMultiArcParametersToEstimate_;

    //! Estimated parameter set with extended multi-arc dynamical parameters only
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > multiArcParametersToEstimate_ ;

    std::shared_ptr< HybridArcResults > variationalPropagationResults_;

    //! Times at which arcs for multi-arc solution start
    std::vector< double > arcStartTimes_;

    //! Function that retrieves the single-arc bodies' initial states as a function of time
    /*!
     *  Function that retrieves the single-arc bodies' initial states as a function of time, is used to update the multi-arc
     *  initial states after the single-arc propagation
     */
    std::function< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >( const double ) >initialStatesFromSingleArcPropagation_;

    double singleArcInitialTime_;


};

//extern template class SingleArcVariationalEquationsSolver< double, double >;
//extern template class MultiArcVariationalEquationsSolver< double, double >;
//extern template class HybridArcVariationalEquationsSolver< double, double >;

} // namespace propagators

} // namespace tudat




#endif // TUDAT_VARIATIONALEQUATIONSSOLVER_H
