/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#include <map>


#include <functional>

#include <Eigen/Core>

#include "tudat/astro/propagators/variationalEquations.h"
#include "tudat/astro/propagators/rotationalMotionQuaternionsStateDerivative.h"

#include "tudat/astro/orbit_determination/acceleration_partials/accelerationPartial.h"


namespace tudat
{

namespace propagators
{

template< typename StateScalarType >
void VariationalEquations::getBodyInitialStatePartialMatrix(
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& stateTransitionAndSensitivityMatrices,
        Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > currentMatrixDerivative )
{
    setBodyStatePartialMatrix( );

    // Add partials of body positions and velocities.
    currentMatrixDerivative.block( 0, 0, totalDynamicalStateSize_, numberOfParameterValues_ ) =
            ( variationalMatrix_.template cast< StateScalarType >( ) * stateTransitionAndSensitivityMatrices );

    if( couplingEntriesToSuppress_ > 0 )
    {
        int numberOfStaticParameters = numberOfParameterValues_ - totalDynamicalStateSize_;
        int numberOfUncoupledEntries = totalDynamicalStateSize_ - couplingEntriesToSuppress_;

        currentMatrixDerivative.block( couplingEntriesToSuppress_, totalDynamicalStateSize_, numberOfUncoupledEntries, numberOfStaticParameters ) =
                variationalMatrix_.template cast< StateScalarType >( ).block(
                    couplingEntriesToSuppress_, couplingEntriesToSuppress_,
                    numberOfUncoupledEntries, numberOfUncoupledEntries ) *
                stateTransitionAndSensitivityMatrices.block(
                    couplingEntriesToSuppress_, totalDynamicalStateSize_, numberOfUncoupledEntries, numberOfStaticParameters );
    }
}

//! Calculates matrix containing partial derivatives of state derivatives w.r.t. body state.
void VariationalEquations::setBodyStatePartialMatrix( )
{
    // Initialize partial matrix
    variationalMatrix_.setZero( );

    if( dynamicalStatesToEstimate_.count( propagators::translational_state ) > 0 )
    {
        int startIndex = stateTypeStartIndices_.at( propagators::translational_state );
        for( unsigned int i = 0; i < dynamicalStatesToEstimate_.at( propagators::translational_state ).size( ); i++ )
        {
            variationalMatrix_.block( startIndex + i * 6, startIndex + i * 6 + 3, 3, 3 ).setIdentity( );
        }
    }

    if( dynamicalStatesToEstimate_.count( propagators::rotational_state ) > 0 )
    {
        Eigen::VectorXd rotationalStates = currentStatesPerTypeInConventionalRepresentation_.at(
                    propagators::rotational_state );

        int startIndex = stateTypeStartIndices_.at( propagators::rotational_state );
        for( unsigned int i = 0; i < dynamicalStatesToEstimate_.at( propagators::rotational_state ).size( ); i++ )
        {
            variationalMatrix_.block( startIndex + i * 7, startIndex + i * 7 , 4, 4 ) =
                    getQuaterionToQuaternionRateMatrix( rotationalStates.segment( 7 * i + 4, 3 ) );
            variationalMatrix_.block( startIndex + i * 7, startIndex + i * 7 + 4, 4, 3 ) =
                    getAngularVelocityToQuaternionRateMatrix( rotationalStates.segment( 7 * i, 4 ) );
        }
    }

    // Iterate over all bodies undergoing accelerations for which initial condition is to be estimated.
    for( std::map< IntegratedStateType, std::vector< std::multimap< std::pair< int, int >,
         std::function< void( Eigen::Block< Eigen::MatrixXd > ) > > > >::iterator
         typeIterator = statePartialList_.begin( ); typeIterator != statePartialList_.end( ); typeIterator++ )
    {
        int startIndex = stateTypeStartIndices_.at( typeIterator->first );
        int currentStateSize = getSingleIntegrationSize( typeIterator->first );
        int entriesToSkipPerEntry = currentStateSize - getGeneralizedAccelerationSize( typeIterator->first );

        for( unsigned int i = 0; i < typeIterator->second.size( ); i++ )
        {
            // Iterate over all bodies exerting an acceleration on this body.
            for( statePartialIterator_ = typeIterator->second.at( i ).begin( );
                 statePartialIterator_ != typeIterator->second.at( i ).end( );
                 statePartialIterator_++ )
            {
                statePartialIterator_->second(
                            variationalMatrix_.block(
                                startIndex + entriesToSkipPerEntry + i* currentStateSize, statePartialIterator_->first.first,
                                currentStateSize - entriesToSkipPerEntry, statePartialIterator_->first.second ) );

            }
        }
    }

    for( unsigned int i = 0; i < statePartialAdditionIndices_.size( ); i++ )
    {
        variationalMatrix_.block( 0, statePartialAdditionIndices_.at( i ).second, totalDynamicalStateSize_, 3 ) +=
                variationalMatrix_.block( 0, statePartialAdditionIndices_.at( i ).first, totalDynamicalStateSize_, 3 );
    }

    for( unsigned int i = 0; i < inertiaTensorsForMultiplication_.size( ); i++ )
    {
        variationalMatrix_.block( inertiaTensorsForMultiplication_.at( i ).first, 0, 3, totalDynamicalStateSize_ ) =
                ( inertiaTensorsForMultiplication_.at( i ).second( ).inverse( ) ) *
                variationalMatrix_.block( inertiaTensorsForMultiplication_.at( i ).first, 0, 3, totalDynamicalStateSize_ ).eval( );
    }

}

//! Function to clear reference/cached values of state derivative partials.
void VariationalEquations::clearPartials( )
{
    for( stateDerivativeTypeIterator_ = stateDerivativePartialList_.begin( );
         stateDerivativeTypeIterator_ != stateDerivativePartialList_.end( );
         stateDerivativeTypeIterator_++ )
    {
        for( unsigned int i = 0; i < stateDerivativeTypeIterator_->second.size( ); i++ )
        {
            for( unsigned int j = 0; j < stateDerivativeTypeIterator_->second.at( i ).size( ); j++ )
            {
                stateDerivativeTypeIterator_->second.at( i ).at( j )->resetCurrentTime( );
            }

        }
    }
}

//! Function (called by constructor) to set up the statePartialList_ member from the state derivative partials
void VariationalEquations::setStatePartialFunctionList( )
{
    std::pair< std::function< void( Eigen::Block< Eigen::MatrixXd > ) >, int > currentDerivativeFunction;


    // Iterate over all state types
    for( std::map< propagators::IntegratedStateType,
         orbit_determination::StateDerivativePartialsMap >::iterator
         stateDerivativeTypeIterator_ = stateDerivativePartialList_.begin( );
         stateDerivativeTypeIterator_ != stateDerivativePartialList_.end( );
         stateDerivativeTypeIterator_++ )
    {
        // Iterate over all bodies undergoing 'accelerations' for which initial state is to be estimated.
        for( unsigned int i = 0; i < stateDerivativeTypeIterator_->second.size( ); i++ )
        {
            std::multimap< std::pair< int, int >, std::function< void( Eigen::Block< Eigen::MatrixXd > ) > >
                    currentBodyPartialList;

            // Iterate over all 'accelerations' from single body on other single body
            for( unsigned int j = 0; j < stateDerivativeTypeIterator_->second.at( i ).size( ); j++ )
            {
                for( std::map< propagators::IntegratedStateType,
                     std::vector< std::pair< std::string, std::string > > >::iterator
                     estimatedStateIterator = dynamicalStatesToEstimate_.begin( );
                     estimatedStateIterator != dynamicalStatesToEstimate_.end( );
                     estimatedStateIterator++ )
                {
                    // Iterate over all bodies to see if body exerting acceleration is also to be estimated (cross-terms)
                    for( unsigned int k = 0; k < estimatedStateIterator->second.size( ); k++ )
                    {
                        currentDerivativeFunction = stateDerivativeTypeIterator_->second.at( i ).at( j )->
                                getDerivativeFunctionWrtStateOfIntegratedBody(
                                    estimatedStateIterator->second.at( k ), estimatedStateIterator->first );

                        // If function is not-empty: add to list.
                        if( currentDerivativeFunction.second != 0 )
                        {
                            currentBodyPartialList.insert(
                                        std::make_pair(
                                            std::make_pair( k * getSingleIntegrationSize( estimatedStateIterator->first ) +
                                                            stateTypeStartIndices_.at( estimatedStateIterator->first ),
                                                            getSingleIntegrationSize( estimatedStateIterator->first ) ),
                                            currentDerivativeFunction.first ) );
                        }
                    }
                }
            }
            statePartialList_[ stateDerivativeTypeIterator_->first ].push_back( currentBodyPartialList );
        }
    }
}


template void VariationalEquations::getBodyInitialStatePartialMatrix< double >(
const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >& stateTransitionAndSensitivityMatrices,
Eigen::Block< Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > > currentMatrixDerivative );

//#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
template void VariationalEquations::getBodyInitialStatePartialMatrix< long double >(
const Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic >& stateTransitionAndSensitivityMatrices,
Eigen::Block< Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic > > currentMatrixDerivative );
//#endif

} // namespace propagators

} // namespace tudat
