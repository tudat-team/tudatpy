/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */



#include "tudat/astro/propagators/stateTransitionMatrixInterface.h"

namespace tudat
{

namespace propagators
{

//! Function to reset the state transition and sensitivity matrix interpolators
void SingleArcCombinedStateTransitionAndSensitivityMatrixInterface::updateMatrixInterpolators(
        const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
        stateTransitionMatrixInterpolator,
        const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
        sensitivityMatrixInterpolator,
        const std::vector< std::pair< int, int > >& statePartialAdditionIndices )
{
    stateTransitionMatrixInterpolator_ = stateTransitionMatrixInterpolator;
    sensitivityMatrixInterpolator_ = sensitivityMatrixInterpolator;

    // Re-order state partial addition indices to match ephemeris update order (inverted in variational equations object)
    statePartialAdditionIndices_.clear( );
    for ( int i = statePartialAdditionIndices.size( ) - 1; i >= 0 ; i-- )
    {
        statePartialAdditionIndices_.push_back( statePartialAdditionIndices[ i ] );
    }
}

//! Function to get the concatenated state transition and sensitivity matrix at a given time.
Eigen::MatrixXd SingleArcCombinedStateTransitionAndSensitivityMatrixInterface::getCombinedStateTransitionAndSensitivityMatrix(
        const double evaluationTime,
        const bool addCentralBodyDependency,
        const std::vector< std::string >& arcDefiningBodies )
{
    combinedStateTransitionMatrix_.setZero( );


    // Set Phi and S matrices.
    combinedStateTransitionMatrix_.block( 0, 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_ ) =
            stateTransitionMatrixInterpolator_->interpolate( evaluationTime );

    if( sensitivityMatrixSize_ > 0 )
    {
        combinedStateTransitionMatrix_.block( 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_, sensitivityMatrixSize_ ) =
                sensitivityMatrixInterpolator_->interpolate( evaluationTime );
    }

    if ( addCentralBodyDependency )
    {
        for( unsigned int i = 0; i < statePartialAdditionIndices_.size( ); i++ )
        {
            combinedStateTransitionMatrix_.block(
                    statePartialAdditionIndices_.at( i ).first, 0, 6, stateTransitionMatrixSize_ + sensitivityMatrixSize_ ) +=
                    combinedStateTransitionMatrix_.block(
                            statePartialAdditionIndices_.at( i ).second, 0, 6, stateTransitionMatrixSize_ + sensitivityMatrixSize_ );
        }
    }

    return combinedStateTransitionMatrix_;
}

}

}

