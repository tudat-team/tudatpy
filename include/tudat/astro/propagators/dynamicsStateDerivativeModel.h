/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DYNAMICSSTATEDERIVATIVEMODEL_H
#define TUDAT_DYNAMICSSTATEDERIVATIVEMODEL_H

#include <map>
#include <utility>

#include <functional>
#include <memory>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/torqueModelTypes.h"
#include "tudat/astro/propagators/bodyMassStateDerivative.h"
#include "tudat/astro/propagators/singleStateTypeDerivative.h"
#include "tudat/astro/propagators/nBodyStateDerivative.h"
#include "tudat/astro/propagators/rotationalMotionStateDerivative.h"
#include "tudat/astro/propagators/variationalEquations.h"

namespace tudat
{

namespace propagators
{

//! Top-level class responsible for single complete function evaluation of dynamics state derivative.
/*!
 *  Top-level class responsible for single complete function evaluation of dynamics state
 *  derivative. This class contains both the EnvironmentUpdater and list of
 *  SingleStateTypeDerivative derived classes that define the full state derivative function, which
 *  fully evaluated by calling the computeStateDerivative function.
 */
template< typename TimeType = double, typename StateScalarType = double >
class DynamicsStateDerivativeModel
{
public:

    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > StateType;
    typedef std::map< IntegratedStateType, std::vector< std::shared_ptr
    < SingleStateTypeDerivative< StateScalarType, TimeType > > > > StateDerivativeCalculatorList;

    //! Derivative model constructor.
    /*!
     *  Derivative model constructor. Takes state derivative model and environment
     *  updater. Constructor checks whether all models use the same environment updater.
     *  \param stateDerivativeModels Vector of state derivative models, with one entry for each type of dynamical equation.
     *  \param environmentUpdateFunction Function which is used to update time-dependent environment models to current time
     *  and state, must be consistent with member environment updaters of stateDerivativeModels entries.
     *  \param variationalEquations Object used for computing the state derivative in the variational equations
     */
    DynamicsStateDerivativeModel(
            const std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > >
            stateDerivativeModels,
            const std::function< void(
                const TimeType, const std::unordered_map< IntegratedStateType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >&,
                const std::vector< IntegratedStateType > ) > environmentUpdateFunction,
            const std::shared_ptr< VariationalEquations > variationalEquations =
            std::shared_ptr< VariationalEquations >( ) ):
        environmentUpdateFunction_( environmentUpdateFunction ), variationalEquations_( variationalEquations ),
        functionEvaluationCounter_( 0 )
    {
        std::vector< IntegratedStateType > stateTypeList;
        totalConventionalStateSize_ = 0;
        totalPropagatedStateSize_ = 0;

        // Iterate over vector of state derivative models, check validity, set member variable map
        // stateDerivativeModels_ and size indices.
        for( unsigned int i = 0; i < stateDerivativeModels.size( ); i++ )
        {
            if( i > 0 )
            {
                if( ( std::find( stateTypeList.begin( ), stateTypeList.end( ),
                                 stateDerivativeModels.at( i )->getIntegratedStateType( ) )
                      != stateTypeList.end( ) )
                        && ( stateDerivativeModels.at( i )->getIntegratedStateType( )
                             != stateDerivativeModels.at( i - 1 )->getIntegratedStateType( ) ) )
                {
                    throw std::runtime_error( "Warning when making hybrid state derivative models, state type " +
                                              std::to_string( stateDerivativeModels.at( i )->getIntegratedStateType( ) )
                                              + " entries are non-contiguous" );
                }
            }

            // Check uniqueness of state derivative type calculator in list
            if( std::find( stateTypeList.begin( ), stateTypeList.end( ),
                           stateDerivativeModels.at( i )->getIntegratedStateType( ) ) == stateTypeList.end( ) )
            {
                conventionalStateTypeSize_[ stateDerivativeModels.at( i )->getIntegratedStateType( ) ] = 0;
                conventionalStateTypeStartIndex_[ stateDerivativeModels.at( i )->getIntegratedStateType( ) ]
                        = totalConventionalStateSize_;

                propagatedStateTypeSize_[ stateDerivativeModels.at( i )->getIntegratedStateType( ) ] = 0;
                propagatedStateTypeStartIndex_[ stateDerivativeModels.at( i )->getIntegratedStateType( ) ]
                        = totalPropagatedStateSize_;
            }
            stateTypeList.push_back( stateDerivativeModels.at( i )->getIntegratedStateType( ) );

            // Set state part sizes
            conventionalStateIndices_[ stateDerivativeModels.at( i )->getIntegratedStateType( ) ].push_back(
                        std::make_pair( totalConventionalStateSize_, stateDerivativeModels.at( i )->getConventionalStateSize( ) ) );
            totalConventionalStateSize_ += stateDerivativeModels.at( i )->getConventionalStateSize( );

            propagatedStateIndices_[ stateDerivativeModels.at( i )->getIntegratedStateType( ) ].push_back(
                        std::make_pair( totalPropagatedStateSize_, stateDerivativeModels.at( i )->getPropagatedStateSize( ) ) );
            totalPropagatedStateSize_ += stateDerivativeModels.at( i )->getPropagatedStateSize( );

            conventionalStateTypeSize_[ stateDerivativeModels.at( i )->getIntegratedStateType( ) ] +=
                    stateDerivativeModels.at( i )->getConventionalStateSize( );
            propagatedStateTypeSize_[ stateDerivativeModels.at( i )->getIntegratedStateType( ) ] +=
                    stateDerivativeModels.at( i )->getPropagatedStateSize( );

            // Set current model in member map.
            stateDerivativeModels_[ stateDerivativeModels.at( i )->getIntegratedStateType( ) ].push_back(
                        stateDerivativeModels.at( i ) );

            currentStatesPerTypeInConventionalRepresentation_[ stateDerivativeModels.at( i )->getIntegratedStateType( )  ] =
                    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero(
                        conventionalStateTypeSize_.at( stateDerivativeModels.at( i )->getIntegratedStateType( )  ), 1 );
        }
    }


    //! Function to calculate the system state derivative
    /*!
     *  Function to calculate the system state derivative, with settings as by last call to
     *  setPropagationSettings function.  Dimensions of state must be consistent with these
     *  settings. Depending on the settings, this function may calculate the dynamical equations
     *  and/or variational equations for a subset of the dynamical equation types that are set in
     *  the stateDerivativeModels_ map.
     *  \param time Current time.
     *  \param state Current complete state.
     *  \return Calculated state derivative.
     */
    StateType computeStateDerivative( const TimeType time, const StateType& state )
    {
        if( !( time == time ) )
        {
            throw std::invalid_argument( "Error when computing system state derivative. Input time is NaN" );
        }

        if( state.hasNaN( ) )
        {
            throw std::invalid_argument( "Error when computing system state derivative. State vector contains NaN" );
        }

        if( !state.allFinite( ) )
        {
            throw std::invalid_argument( "Error when computing system state derivative. State vector contains Inf" );
        }
//        std::cout << "Computing state derivative: " <<time<<" "<<state.transpose( ) << std::endl;

        // Initialize state derivative
        if( stateDerivative_.rows( ) != state.rows( ) || stateDerivative_.cols( ) != state.cols( )  )
        {
            stateDerivative_.resize( state.rows( ), state.cols( ) );
        }

        // If dynamical equations are integrated, update the environment with the current state.
        if( evaluateDynamicsEquations_ )
        {
            // Iterate over all types of equations.
            for( stateDerivativeModelsIterator_ = stateDerivativeModels_.begin( );
                 stateDerivativeModelsIterator_ != stateDerivativeModels_.end( );
                 stateDerivativeModelsIterator_++ )
            {
                for( unsigned int i = 0; i < stateDerivativeModelsIterator_->second.size( ); i++ )
                {
                    stateDerivativeModelsIterator_->second.at( i )->clearStateDerivativeModel( );
                }
            }

            convertCurrentStateToGlobalRepresentationPerType( state, time, evaluateVariationalEquations_ );
            environmentUpdateFunction_( time, currentStatesPerTypeInConventionalRepresentation_,
                                        integratedStatesFromEnvironment_ );
        }
        else
        {
            environmentUpdateFunction_(
                        time, std::unordered_map<
                        IntegratedStateType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >( ),
                        integratedStatesFromEnvironment_ );
        }

        if( evaluateVariationalEquations_ )
        {
            variationalEquations_->clearPartials( );
        }

        // If dynamical equations are integrated, evaluate dynamics state derivatives.
        std::pair< int, int > currentIndices;
        if( evaluateDynamicsEquations_ )
        {
            // Iterate over all types of equations.
            for( stateDerivativeModelsIterator_ = stateDerivativeModels_.begin( );
                 stateDerivativeModelsIterator_ != stateDerivativeModels_.end( );
                 stateDerivativeModelsIterator_++ )
            {
                for( unsigned int i = 0; i < stateDerivativeModelsIterator_->second.size( ); i++ )
                {
                    // Update state derivative models
                    stateDerivativeModelsIterator_->second.at( i )->updateStateDerivativeModel( time );
                }
            }

            for( stateDerivativeModelsIterator_ = stateDerivativeModels_.begin( );
                 stateDerivativeModelsIterator_ != stateDerivativeModels_.end( );
                 stateDerivativeModelsIterator_++ )
            {
                for( unsigned int i = 0; i < stateDerivativeModelsIterator_->second.size( ); i++ )
                {
                    // Evaluate and set current dynamical state derivative
                    currentIndices = propagatedStateIndices_.at( stateDerivativeModelsIterator_->first ).at( i );

                    stateDerivativeModelsIterator_->second.at( i )->calculateSystemStateDerivative(
                                time, state.block( currentIndices.first, dynamicsStartColumn_, currentIndices.second, 1 ),
                                stateDerivative_.block( currentIndices.first, dynamicsStartColumn_, currentIndices.second, 1 ) );
                }
            }
        }

        // If variational equations are to be integrated: evaluate and set.
        if( evaluateVariationalEquations_ )
        {
            variationalEquations_->updatePartials( time, currentStatesPerTypeInConventionalRepresentation_ );

            variationalEquations_->evaluateVariationalEquations< StateScalarType >(
                        time, state.block( 0, 0, totalConventionalStateSize_, variationalEquations_->getNumberOfParameterValues( ) ),
                        stateDerivative_.block( 0, 0, totalConventionalStateSize_, variationalEquations_->getNumberOfParameterValues( ) ) );
        }

        // Update counters
        functionEvaluationCounter_++;
        cumulativeFunctionEvaluationCounter_[ time ] = functionEvaluationCounter_;

        return stateDerivative_;

    }

    //! Function to calculate the system state derivative with double precision, regardless of template arguments.
    /*!
     *  Function to calculate the system state derivative with double precision, regardless of template arguments.
     *  \sa computeStateDerivative
     *  \param time Current time.
     *  \param state Current complete state.
     *  \return Calculated state derivative.
     */
    Eigen::MatrixXd computeStateDoubleDerivative(
            const double time, const Eigen::MatrixXd& state )
    {
        return computeStateDerivative( static_cast< TimeType >( time ),
                                       state.template cast< StateScalarType >( ) ).template cast< double >( );
    }

    //! Function to convert the state in the conventional form to the propagator-specific form.
    /*!
     * Function to convert the state in the conventional form to the propagator-specific form.  The
     * conventional form is one that is typically used to represent the current state in the
     * environment (e.g. Body class). For translational dynamics this is the Cartesian position and
     * velocity).
     * \param outputState State in 'conventional form'
     * \param time Current time at which the state is valid.
     * \return State (outputState), converted to the 'propagator-specific form'
     */
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& outputState,
            const TimeType& time )
    {
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > internalState =
                Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( totalPropagatedStateSize_, 1 );

        // Iterate over all state derivative models and convert associated state entries
        std::vector< std::pair< int, int > > currentConventionalStateIndices;
        std::vector< std::pair< int, int > > currentPropagatedStateIndices;
        for( stateDerivativeModelsIterator_ = stateDerivativeModels_.begin( );
             stateDerivativeModelsIterator_ != stateDerivativeModels_.end( );
             stateDerivativeModelsIterator_++ )
        {
            currentConventionalStateIndices = conventionalStateIndices_.at( stateDerivativeModelsIterator_->first );
            currentPropagatedStateIndices = propagatedStateIndices_.at( stateDerivativeModelsIterator_->first );
            for( unsigned int i = 0; i < stateDerivativeModelsIterator_->second.size( ); i++ )
            {
                internalState.segment( currentPropagatedStateIndices.at( i ).first,
                                       currentPropagatedStateIndices.at( i ).second ) =
                        stateDerivativeModelsIterator_->second.at( i )->convertFromOutputSolution(
                            outputState.segment( currentConventionalStateIndices.at( i ).first,
                                                 currentConventionalStateIndices.at( i ).second ), time );
            }
        }

        return internalState;
    }

    //! Function to convert the propagator-specific form of the state to the conventional form.
    /*!
     * Function to convert the propagator-specific form of the state to the conventional form. The
     * conventional form is one that is typically used to represent the current state in the
     * environment (e.g. Body class). For translational dynamics this is the Cartesian position and
     * velocity).  In contrast to the convertCurrentStateToGlobalRepresentation function, this
     * function does not provide the state in the inertial frame, but instead provides it in the
     * frame in which it is propagated.  \param internalSolution State in propagator-specific form
     * (i.e. form that is used in numerical integration).
     * \param time Current time at which the state is valid.
     * \return State (internalSolution), converted to the 'conventional form'
     */
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > convertToOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& internalSolution,
            const TimeType& time )
    {
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > outputState =
                Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( totalConventionalStateSize_, 1 );

        // Iterate over all state derivative models and convert associated state entries
        std::vector< std::pair< int, int > > currentConventionalStateIndices;
        std::vector< std::pair< int, int > > currentPropagatedStateIndices;
        for( stateDerivativeModelsIterator_ = stateDerivativeModels_.begin( );
             stateDerivativeModelsIterator_ != stateDerivativeModels_.end( );
             stateDerivativeModelsIterator_++ )
        {
            currentConventionalStateIndices = conventionalStateIndices_.at( stateDerivativeModelsIterator_->first );
            currentPropagatedStateIndices = propagatedStateIndices_.at( stateDerivativeModelsIterator_->first );
            for( unsigned int i = 0; i < stateDerivativeModelsIterator_->second.size( ); i++ )
            {
                stateDerivativeModelsIterator_->second.at( i )->convertToOutputSolution(
                            internalSolution.segment(
                                currentPropagatedStateIndices.at( i ).first,
                                currentPropagatedStateIndices.at( i ).second ), time,
                            outputState.block( currentConventionalStateIndices.at( i ).first, 0,
                                               currentConventionalStateIndices.at( i ).second, 1 ) );
            }
        }

        return outputState;
    }

    //! Function to convert a state history from propagator-specific form to the conventional form.
    /*!
     * Function to convert a state history from propagator-specific form to the conventional form
     * (not necessarily in inertial frame).
     * \sa DynamicsStateDerivativeModel::convertToOutputSolution
     * \param convertedSolution State history (rawSolution), converted to the 'conventional form' (by reference)
     * \param rawSolution State history in propagator-specific form (i.e. form that is used in
     *        numerical integration).
     */
    void convertNumericalStateSolutionsToOutputSolutions(
            std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& convertedSolution,
            const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& rawSolution )
    {
        // Iterate over all times.
        for( typename std::map< TimeType, Eigen::Matrix< StateScalarType,
             Eigen::Dynamic, 1 > >::const_iterator
             stateIterator = rawSolution.begin( ); stateIterator != rawSolution.end( ); stateIterator++ )
        {
            // Convert solution at this time to output (Cartesian with propagation origin frame for
            // translational dynamics) solution
            convertedSolution[ stateIterator->first ] =
                    convertToOutputSolution( stateIterator->second, stateIterator->first );
        }
    }

    //! Function to process the state vector during propagation.
    /*!
     * Function to process the state vector during propagation.
     * \param unprocessedState State before processing.
     * \return Processed state (returned by reference).
     */
    void postProcessState( Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& unprocessedState )
    {
        // Iterate over all state derivative models and post-process associated state entries
        std::vector< std::pair< int, int > > currentIndices;
        for( stateDerivativeModelsIterator_ = stateDerivativeModels_.begin( );
             stateDerivativeModelsIterator_ != stateDerivativeModels_.end( );
             stateDerivativeModelsIterator_++ )
        {
            currentIndices = propagatedStateIndices_.at( stateDerivativeModelsIterator_->first );
            for( unsigned int i = 0; i < stateDerivativeModelsIterator_->second.size( ); i++ )
            {
                if ( stateDerivativeModelsIterator_->second.at( i )->isStateToBePostProcessed( ) )
                {
                    stateDerivativeModelsIterator_->second.at( i )->postProcessState(
                                unprocessedState.block( currentIndices.at( i ).first, 0,
                                                        currentIndices.at( i ).second, 1 ) );
                }
            }
        }
    }

    //! Function to process the state vector and variational equations during propagation.
    /*!
     * Function to process the state vector and variational equations during propagation.
     * \param unprocessedState State before processing.
     * \return Processed state (returned by reference).
     */
    void postProcessStateAndVariationalEquations(
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& unprocessedState )
    {
        // Iterate over all state derivative models and post-process associated state entries
        std::vector< std::pair< int, int > > currentIndices;
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > currentUnprocessedState;
        for( stateDerivativeModelsIterator_ = stateDerivativeModels_.begin( );
             stateDerivativeModelsIterator_ != stateDerivativeModels_.end( );
             stateDerivativeModelsIterator_++ )
        {
            currentIndices = propagatedStateIndices_.at( stateDerivativeModelsIterator_->first );
            for( unsigned int i = 0; i < stateDerivativeModelsIterator_->second.size( ); i++ )
            {
            	if ( stateDerivativeModelsIterator_->second.at( i )->isStateToBePostProcessed( ) )
                {
                    currentUnprocessedState = unprocessedState.block( currentIndices.at( i ).first, dynamicsStartColumn_,
                                                                      currentIndices.at( i ).second, 1 );
                    stateDerivativeModelsIterator_->second.at( i )->postProcessState( currentUnprocessedState );
                    unprocessedState.block( currentIndices.at( i ).first, dynamicsStartColumn_,
                                            currentIndices.at( i ).second, 1 ) = currentUnprocessedState;
               	}
            }
        }
    }

    //! Function to add variational equations to the state derivative model
    /*!
     * Function to add variational equations to the state derivative model.
     * \param variationalEquations Object used for computing the state derivative in the variational equations
     */
    void addVariationalEquations( std::shared_ptr< VariationalEquations > variationalEquations )
    {
        variationalEquations_ = variationalEquations;
    }

    //! Function to set which segments of the full state to propagate
    /*!
     * Function to set which segments of the full state to propagate, i.e. whether to propagate the
     * variational/dynamical equations, and which types of the dynamics to propagate.
     * \param stateTypesToNotIntegrate Types of dynamics to propagate
     * \param evaluateDynamicsEquations Boolean to denote whether the dynamical equations are to be propagated or not
     * \param evaluateVariationalEquations Boolean to denote whether the variational equations are to be propagated or not
     */
    void setPropagationSettings(
            const std::vector< IntegratedStateType >& stateTypesToNotIntegrate,
            const bool evaluateDynamicsEquations,
            const bool evaluateVariationalEquations )
    {
        integratedStatesFromEnvironment_ = stateTypesToNotIntegrate;
        evaluateDynamicsEquations_ = evaluateDynamicsEquations;
        evaluateVariationalEquations_ = evaluateVariationalEquations;

        if( evaluateVariationalEquations_ )
        {
            dynamicsStartColumn_ = variationalEquations_->getNumberOfParameterValues( );
        }
        else
        {
            dynamicsStartColumn_ = 0;
        }
    }

    //! Function to update the settings of the state derivative models with new initial states
    /*!
     * Function to update the settings of the state derivative models with new initial states. This function is
     * called when using, for instance and Encke propagator for the translational dynamics, and the reference orbits
     * are modified.
     * \param initialBodyStates New initial state for the full propagated dynamics.
     */
    void updateStateDerivativeModelSettings(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialBodyStates )
    {
        // Iterate over all dynamics types
        for( stateDerivativeModelsIterator_ = stateDerivativeModels_.begin( );
             stateDerivativeModelsIterator_ != stateDerivativeModels_.end( );
             stateDerivativeModelsIterator_++ )
        {
            switch( stateDerivativeModelsIterator_->first )
            {
            case translational_state:
            {
                for( unsigned int i = 0; i < stateDerivativeModelsIterator_->second.size( ); i++ )
                {
                    std::shared_ptr< NBodyStateDerivative< StateScalarType, TimeType > > currentTranslationalStateDerivative =
                            std::dynamic_pointer_cast< NBodyStateDerivative< StateScalarType, TimeType > >(
                                stateDerivativeModelsIterator_->second.at( i ) );
                    switch( currentTranslationalStateDerivative->getTranslationalPropagatorType( ) )
                    {
                    case cowell:
                        break;
                    case encke:
                        throw std::runtime_error( "Error, reference orbit not reset in Encke propagator" );
                        break;
                    case gauss_keplerian:
                        break;
                    case gauss_modified_equinoctial:
                        break;
                    case unified_state_model_quaternions:
                        break;
                    case unified_state_model_modified_rodrigues_parameters:
                        break;
                    case unified_state_model_exponential_map:
                        break;
                    default:
                        throw std::runtime_error( "Error when updating state derivative model settings, did not recognize translational propagator type" );
                        break;
                    }
                }
            }
            case rotational_state:
                break;
            case body_mass_state:
                break;
            case custom_state:
                break;
            default:
                throw std::runtime_error( "Error when updating state derivative model settings, did not recognize dynamics type" );
                break;
            }
        }
    }

    //! Function to get complete list of state derivative models, sorted per state type.
    /*!
     * Function to get complete list of state derivative models, sorted per state type.
     * \return Complete list of state derivative models, sorted per state type.
     */
    std::unordered_map< IntegratedStateType, std::vector< std::shared_ptr
    < SingleStateTypeDerivative< StateScalarType, TimeType > > > > getStateDerivativeModels( )
    {
        return stateDerivativeModels_;
    }

    //! Function to get state start index per state type in the complete state vector.
    /*!
     * Function to get state start index per state type in the complete state vector.
     * \return State start index per state type in the complete state vector.
     */
    std::map< IntegratedStateType, int > getStateTypeStartIndices( )
    {
        return conventionalStateTypeStartIndex_;
    }

    //! Function to retrieve number of calls to the computeStateDerivative function
    /*!
     * Function to retrieve number of calls to the computeStateDerivative function since object creation/last call to
     * resetFunctionEvaluationCounter function
     * \return Number of calls to the computeStateDerivative function since object creation/last call to
     * resetFunctionEvaluationCounter function
     */
    unsigned int getNumberOfFunctionEvaluations( )
    {
        return functionEvaluationCounter_;
    }

    //! Function to reset the number of calls to the computeStateDerivative function to zero.
    /*!
     * Function to resetr the number of calls to the computeStateDerivative function to zero.  Typically called before any
     * start of numerical integration of dynamics (automatically by DynamicsSimulator)
     */
    void resetFunctionEvaluationCounter( )
    {
        functionEvaluationCounter_ = 0;
    }

    //! Function to retrieve number of calls to the computeStateDerivative function per time step
    /*!
     * Function to retrieve number of calls to the computeStateDerivative function per time step since object
     * reation/last call to resetFunctionEvaluationCounter function
     * \return Number of calls to the computeStateDerivative function since object creation/last call to
     * resetFunctionEvaluationCounter function
     */
    std::map< TimeType, unsigned int > getCumulativeNumberOfFunctionEvaluations( )
    {
        return cumulativeFunctionEvaluationCounter_;
    }

    //! Function to reset the number of calls to the computeStateDerivative function to zero.
    /*!
     * Function to resetr the number of calls to the computeStateDerivative function to zero.  Typically called before any
     * start of numerical integration of dynamics (automatically by DynamicsSimulator)
     */
    void resetCumulativeFunctionEvaluationCounter( )
    {
        cumulativeFunctionEvaluationCounter_.clear( );
    }

    //! Function to retrieve the object used for computing the state derivative in the variational equations
    /*!
     * Function to retrieve the object used for computing the state derivative in the variational equations
     * \return Object used for computing the state derivative in the variational equations
     */
    std::shared_ptr< VariationalEquations > getVariationalEquationsCalculator( )
    {
        return variationalEquations_;
    }


private:

    //! Function to convert the to the conventional form in the global frame per dynamics type.
    /*!
     * Function to convert the propagator-specific form of the state to the conventional form in the global frame, split
     * by dynamics type. This function updates the currentStatesPerTypeInConventionalRepresentation_ to the current state
     * and time.
     * The conventional form is one that is typically used to represent the current state in the environment
     * (e.g. Body class). For translational dynamics this is the Cartesian position and velocity).
     * The inertial frame is typically the barycenter with J2000/ECLIPJ2000 orientation, but may differ depending on
     * simulation settings
     * \param state State in propagator-specific form (i.e. form that is used in numerical integration).
     * \param time Current time at which the state is valid.
     * \param stateIncludesVariationalState Boolean defining whether the stae includes the state transition/sensitivity
     * matrices
     */
    void convertCurrentStateToGlobalRepresentationPerType(
            const StateType& state, const TimeType& time, const bool stateIncludesVariationalState )
    {
        int startColumn = 0;
        if( stateIncludesVariationalState )
        {
            startColumn = variationalEquations_->getNumberOfParameterValues( );
        }
        else
        {
            startColumn = 0;
        }

        std::pair< int, int > currentPropagatedIndices, currentConventionalIndices;

        // Iterate over all state derivative models
        for( stateDerivativeModelsIterator_ = stateDerivativeModels_.begin( );
             stateDerivativeModelsIterator_ != stateDerivativeModels_.end( );
             stateDerivativeModelsIterator_++ )
        {
            int currentStateTypeSize = 0;

            // Iterate over all state derivative models of current type
            for( unsigned int i = 0; i < stateDerivativeModelsIterator_->second.size( ); i++ )
            {
                // Get state block indices of current state derivative model
                currentPropagatedIndices = propagatedStateIndices_.at( stateDerivativeModelsIterator_->first ).at( i );
                currentConventionalIndices = conventionalStateIndices_.at( stateDerivativeModelsIterator_->first ).at( i );

                // Set current block in split state (in global form)
                stateDerivativeModelsIterator_->second.at( i )->convertCurrentStateToGlobalRepresentation(
                            state.block( currentPropagatedIndices.first, startColumn, currentPropagatedIndices.second, 1 ), time,
                            currentStatesPerTypeInConventionalRepresentation_.at(
                                stateDerivativeModelsIterator_->first ).block(
                                currentStateTypeSize, 0, currentConventionalIndices.second, 1 ) );

            }
        }
    }

    std::function<
    void( const TimeType, const std::unordered_map< IntegratedStateType,
          Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >&,
          const std::vector< IntegratedStateType > ) > environmentUpdateFunction_;

    //! Object used for computing the state derivative in the variational equations
    std::shared_ptr< VariationalEquations > variationalEquations_;

    //! Map that denotes for each state derivative model the start index and size of the associated
    //! state in the full state vector.
    std::map< IntegratedStateType, std::vector< std::pair< int, int > > > conventionalStateIndices_;

    std::map< IntegratedStateType, std::vector< std::pair< int, int > > > propagatedStateIndices_;

    //! State size per state type in the complete state vector.
    std::map< IntegratedStateType, int > conventionalStateTypeSize_;

    std::map< IntegratedStateType, int > propagatedStateTypeSize_;

    //! State start index per state type in the complete state vector.
    std::map< IntegratedStateType, int > conventionalStateTypeStartIndex_;

    std::map< IntegratedStateType, int > propagatedStateTypeStartIndex_;

    //! Complete list of state derivative models, sorted per state type.
    std::unordered_map< IntegratedStateType,
    std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > > stateDerivativeModels_;

    //! Predefined iterator for computational efficiency.
    typename std::unordered_map< IntegratedStateType, std::vector< std::shared_ptr
    < SingleStateTypeDerivative< StateScalarType, TimeType > > > >::iterator stateDerivativeModelsIterator_;

    //! Total length of conventional state vector.
    /*!
     *  Total length of conventional state vector. For instance, for translational propagation, this is the
     *  length of the Cartesian state, whereas for rotational propagation, it is the length of the quaternion
     *  state.
     */
    int totalConventionalStateSize_;

    //! Total length of propagated state vector.
    int totalPropagatedStateSize_;

    //! List of states that are not propagated in current numerical integration, i.e, for which
    //! current state is taken from the environment.
    std::vector< IntegratedStateType > integratedStatesFromEnvironment_;

    //! Boolean denoting whether the equations of motion are to be propagated or not.
    bool evaluateDynamicsEquations_;

    //! Boolean denoting whether the variational equations are to be propagated or not.
    bool evaluateVariationalEquations_;

    //! Start index in propagated matrix of the equations of motion (=0 if variational equations are
    //! not propagated).
    int dynamicsStartColumn_;

    //! Current state derivative, as computed by computeStateDerivative.
    StateType stateDerivative_;

    //! Current state in 'conventional' representation, computed from current propagated state by
    //! convertCurrentStateToGlobalRepresentationPerType
    std::unordered_map< IntegratedStateType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >
    currentStatesPerTypeInConventionalRepresentation_;

    //! Variable to keep track of the number of calls to the computeStateDerivative function
    unsigned int functionEvaluationCounter_ = 0;

    //! Variable to keep track of the number of calls to the computeStateDerivative function per time step
    std::map< TimeType, unsigned int > cumulativeFunctionEvaluationCounter_;
};

extern template class DynamicsStateDerivativeModel< double, double >;

#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
extern template class DynamicsStateDerivativeModel< Time, double >;
extern template class DynamicsStateDerivativeModel< double, long double >;
extern template class DynamicsStateDerivativeModel< Time, long double >;
#endif


//! Function to retrieve a single given acceleration model from a list of models
/*!
 *  Function to retrieve a single given acceleration model, determined by
 *  the body exerting and undergoing the acceleration, as well as the acceleration type, from a list of
 *  state derivative models.
 *  \param bodyUndergoingAcceleration Name of body undergoing the acceleration.
 *  \param bodyExertingAcceleration Name of body exerting the acceleration.
 *  \param stateDerivativeModels Complete list of state derivativ models
 *  \param accelerationModelType Type of acceleration model that is to be retrieved.
 */
template< typename TimeType = double, typename StateScalarType = double >
std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel3d > > getAccelerationBetweenBodies(
        const std::string bodyUndergoingAcceleration,
        const std::string bodyExertingAcceleration,
        const std::unordered_map< IntegratedStateType,
        std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >& stateDerivativeModels,
        const basic_astrodynamics::AvailableAcceleration accelerationModelType )
{
    std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
            listOfSuitableAccelerationModels;

    // Retrieve acceleration models
    if( stateDerivativeModels.count( propagators::translational_state ) == 1 )
    {
        basic_astrodynamics::AccelerationMap accelerationModelList =
                std::dynamic_pointer_cast< NBodyStateDerivative< StateScalarType, TimeType > >(
                    stateDerivativeModels.at( propagators::translational_state ).at( 0 ) )->getFullAccelerationsMap( );

        if( accelerationModelList.count( bodyUndergoingAcceleration ) == 0 )
        {

            std::string errorMessage = "Error when getting acceleration between bodies, no translational dynamics models acting on " +
                    bodyUndergoingAcceleration + " are found";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            // Retrieve accelerations acting on bodyUndergoingAcceleration
            if( accelerationModelList.at( bodyUndergoingAcceleration ).count( bodyExertingAcceleration ) == 0 )
            {
                std::string errorMessage = "Error when getting acceleration between bodies, no translational dynamics models by " +
                        bodyExertingAcceleration + " acting on " + bodyUndergoingAcceleration + " are found";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                // Retrieve required acceleration.
                listOfSuitableAccelerationModels = basic_astrodynamics::getAccelerationModelsOfType(
                            accelerationModelList.at( bodyUndergoingAcceleration ).at( bodyExertingAcceleration ), accelerationModelType );
            }
        }
    }
    else
    {
        std::string errorMessage = "Error when getting acceleration between bodies, no translational dynamics models found";
        throw std::runtime_error( errorMessage );
    }
    return listOfSuitableAccelerationModels;
}

//! Function to retrieve a single given torque model from a list of models
/*!
 *  Function to retrieve a single given torque model, determined by
 *  the body exerting and undergoing the torque, as well as the torque type, from a list of
 *  state derivative models.
 *  \param bodyUndergoingTorque Name of body undergoing the torque.
 *  \param bodyExertingTorque Name of body exerting the torque.
 *  \param stateDerivativeModels Complete list of state derivativ models
 *  \param torqueModelType Type of torque model that is to be retrieved.
 */
template< typename TimeType = double, typename StateScalarType = double >
std::vector< std::shared_ptr< basic_astrodynamics::TorqueModel > > getTorqueBetweenBodies(
        const std::string bodyUndergoingTorque,
        const std::string bodyExertingTorque,
        const std::unordered_map< IntegratedStateType,
        std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >& stateDerivativeModels,
        const basic_astrodynamics::AvailableTorque torqueModelType )
{
    std::vector< std::shared_ptr< basic_astrodynamics::TorqueModel > >
            listOfSuitableTorqueModels;

    // Retrieve torque models
    if( stateDerivativeModels.count( propagators::rotational_state ) == 1 )
    {
        basic_astrodynamics::TorqueModelMap torqueModelList =
                std::dynamic_pointer_cast< RotationalMotionStateDerivative< StateScalarType, TimeType > >(
                    stateDerivativeModels.at( propagators::rotational_state ).at( 0 ) )->getTorquesMap( );
        if( torqueModelList.count( bodyUndergoingTorque ) == 0 )
        {
            std::string errorMessage = "Error when getting torque between bodies, no translational dynamics models acting on " +
                    bodyUndergoingTorque + " are found";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            // Retrieve torques acting on bodyUndergoingTorque
            if( torqueModelList.at( bodyUndergoingTorque ).count( bodyExertingTorque ) == 0 )
            {
                std::string errorMessage = "Error when getting torque between bodies, no translational dynamics models by " +
                        bodyExertingTorque + " acting on " + bodyUndergoingTorque + " are found";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                // Retrieve required torque.
                listOfSuitableTorqueModels = basic_astrodynamics::getTorqueModelsOfType(
                            torqueModelList.at( bodyUndergoingTorque ).at( bodyExertingTorque ), torqueModelType );
            }
        }
    }
    else
    {
        std::string errorMessage = "Error when getting torque between bodies, no translational dynamics models found";
        throw std::runtime_error( errorMessage );
    }
    return listOfSuitableTorqueModels;
}

//! Function to retrieve the state derivative models for translational dynamics of given body.
/*!
 *  Function to retrieve the state derivative models for translational dynamics (object of derived class from
 *  NBodyStateDerivative) of given body from full list of state derivative models
 *  \param bodyUndergoingAcceleration Name of body for which state derivative model is to be retrieved
 *  \param stateDerivativeModels Complete list of state derivativ models
 *  \return State derivative model computing derivative of translational state of requested body.
 */
template< typename TimeType = double, typename StateScalarType = double >
std::shared_ptr< NBodyStateDerivative< StateScalarType, TimeType > > getTranslationalStateDerivativeModelForBody(
        const std::string bodyUndergoingAcceleration,
        const std::unordered_map< IntegratedStateType,
        std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >& stateDerivativeModels )

{
    bool modelFound = 0;
    std::shared_ptr< NBodyStateDerivative< StateScalarType, TimeType > > modelForBody;

    // Check if translational state derivative models exists
    if( stateDerivativeModels.count( propagators::translational_state ) > 0 )
    {
        for( unsigned int i = 0; i < stateDerivativeModels.at( propagators::translational_state ).size( ); i++ )
        {
            std::shared_ptr< NBodyStateDerivative< StateScalarType, TimeType > > nBodyModel =
                    std::dynamic_pointer_cast< NBodyStateDerivative< StateScalarType, TimeType > >(
                        stateDerivativeModels.at( propagators::translational_state ).at( i ) );
            std::vector< std::string > propagatedBodies = nBodyModel->getBodiesToBeIntegratedNumerically( );

            // Check if bodyUndergoingAcceleration is propagated by bodyUndergoingAcceleration
            if( std::find( propagatedBodies.begin( ), propagatedBodies.end( ), bodyUndergoingAcceleration )
                    != propagatedBodies.end( ) )
            {
                if( modelFound == true )
                {
                    std::string errorMessage = "Error when getting translational dynamics model for " +
                            bodyUndergoingAcceleration + ", multiple models found";
                    throw std::runtime_error( errorMessage );
                }
                else
                {
                    modelForBody = nBodyModel;
                    modelFound = true;
                }
            }
        }
    }
    else
    {
        std::string errorMessage = "Error when getting translational dynamics model for " +
                bodyUndergoingAcceleration + " no translational dynamics models found";
        throw std::runtime_error( errorMessage );
    }
    return modelForBody;
}

template< typename TimeType = double, typename StateScalarType = double >
std::shared_ptr< RotationalMotionStateDerivative< StateScalarType, TimeType > > getRotationalStateDerivativeModelForBody(
        const std::string bodyUndergoingTorque,
        const std::unordered_map< IntegratedStateType,
        std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >& stateDerivativeModels )

{
    bool modelFound = 0;
    std::shared_ptr< RotationalMotionStateDerivative< StateScalarType, TimeType > > modelForBody;

    // Check if translational state derivative models exists
    if( stateDerivativeModels.count( propagators::rotational_state ) > 0 )
    {
        for( unsigned int i = 0; i < stateDerivativeModels.at( propagators::rotational_state ).size( ); i++ )
        {
            std::shared_ptr< RotationalMotionStateDerivative< StateScalarType, TimeType > > rotationalDynamicsModel =
                    std::dynamic_pointer_cast< RotationalMotionStateDerivative< StateScalarType, TimeType > >(
                        stateDerivativeModels.at( propagators::rotational_state ).at( i ) );
            std::vector< std::string > propagatedBodies = rotationalDynamicsModel->getBodiesToBeIntegratedNumerically( );

            // Check if bodyUndergoingTorque is propagated by bodyUndergoingTorque
            if( std::find( propagatedBodies.begin( ), propagatedBodies.end( ), bodyUndergoingTorque )
                    != propagatedBodies.end( ) )
            {
                if( modelFound == true )
                {
                    std::string errorMessage = "Error when getting rotational dynamics model for " +
                            bodyUndergoingTorque + ", multiple models found";
                    throw std::runtime_error( errorMessage );
                }
                else
                {
                    modelForBody = rotationalDynamicsModel;
                    modelFound = true;
                }
            }
        }
    }
    else
    {
        std::string errorMessage = "Error when getting translational dynamics model for " +
                bodyUndergoingTorque + " no translational dynamics models found";
        throw std::runtime_error( errorMessage );
    }
    return modelForBody;
}

//! Function to retrieve the state derivative models for body mass of given body.
/*!
 *  Function to retrieve the state derivative models for body mass (object of class BodyMassStateDerivative)
 *  of given body from full list of state derivative models
 *  \param bodyWithMassDerivative Name of body for which mass derivative model is to be retrieved
 *  \param stateDerivativeModels Complete list of state derivative models
 *  \return State derivative model computing derivative of body mass of requested body.
 */
template< typename TimeType = double, typename StateScalarType = double >
std::shared_ptr< BodyMassStateDerivative< StateScalarType, TimeType > > getBodyMassStateDerivativeModelForBody(
        const std::string bodyWithMassDerivative,
        const std::unordered_map< IntegratedStateType,
        std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >& stateDerivativeModels )

{
    bool modelFound = 0;
    std::shared_ptr< BodyMassStateDerivative< StateScalarType, TimeType > > modelForBody;

    // Check if body mass derivative models exists
    if( stateDerivativeModels.count( propagators::body_mass_state ) > 0 )
    {
        for( unsigned int i = 0; i < stateDerivativeModels.at( propagators::body_mass_state ).size( ); i++ )
        {
            // Test consistency of current object
            std::shared_ptr< BodyMassStateDerivative< StateScalarType, TimeType > > massRateModel =
                    std::dynamic_pointer_cast< BodyMassStateDerivative< StateScalarType, TimeType > >(
                        stateDerivativeModels.at( propagators::body_mass_state ).at( i ) );
            if( massRateModel == nullptr )
            {
                throw std::runtime_error( "Error when finding body's mass rate model, model type is inconsistent" );
            }

            std::vector< std::string > propagatedBodies = massRateModel->getBodiesToIntegrate( );

            // Check if bodyWithMassDerivative is propagated by massRateModel
            if( std::find( propagatedBodies.begin( ), propagatedBodies.end( ), bodyWithMassDerivative )
                    != propagatedBodies.end( ) )
            {
                if( modelFound == true )
                {
                    std::string errorMessage = "Error when getting mass rate  model for " +
                            bodyWithMassDerivative + ", multiple models found";
                    throw std::runtime_error( errorMessage );
                }
                else
                {
                    modelForBody = massRateModel;
                    modelFound = true;
                }
            }
        }
    }
    else
    {
        std::string errorMessage = "Error when getting mass rate model for " +
                bodyWithMassDerivative + " no mass rate models found";
        throw std::runtime_error( errorMessage );
    }

    return modelForBody;
}

#if(TUDAT_BUILD_WITH_ESTIMATION_TOOLS )
//! Function to retrieve specific acceleration partial object from list of state derivative partials
/*!
 * Function to retrieve specific acceleration partial object from list of state derivative partials
 * \param accelerationPartials List of state derivative partials
 * \param accelerationType Type of acceleration for which partial is to be retrieved
 * \param bodyExertingAcceleration Name of body exerting acceleration
 * \param bodyUndergoignAcceleration Name of body undergoing acceleration
 * \param centralBody Name of central body (only required for third-body acceleration types).
 * \return Required acceleration partial from list.
 */
std::shared_ptr< acceleration_partials::AccelerationPartial > getAccelerationPartialForBody(
        const orbit_determination::StateDerivativePartialsMap& accelerationPartials,
        const basic_astrodynamics::AvailableAcceleration accelerationType,
        const std::string& bodyExertingAcceleration,
        const std::string& bodyUndergoignAcceleration,
        const std::string& centralBody = "" );
#endif

template< typename TimeType = double, typename StateScalarType = double,
          typename ConversionClassType = DynamicsStateDerivativeModel< TimeType, StateScalarType > >
void convertNumericalStateSolutionsToOutputSolutions(
        std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& convertedSolution,
        const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& rawSolution,
        std::shared_ptr< ConversionClassType > converterClass )
{
    // Iterate over all times.
    for( typename std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >::const_iterator stateIterator =
         rawSolution.begin( ); stateIterator != rawSolution.end( ); stateIterator++ )
    {
        // Convert solution at this time to output (typically ephemeris frame of given body) solution
        convertedSolution[ stateIterator->first ] = converterClass->convertToOutputSolution( stateIterator->second, stateIterator->first );
    }
}

} // namespace propagators

} // namespace tudat


#endif // TUDAT_DYNAMICSSTATEDERIVATIVEMODEL_H
