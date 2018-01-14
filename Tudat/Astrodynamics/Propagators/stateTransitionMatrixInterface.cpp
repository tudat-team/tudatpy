/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/Propagators/stateTransitionMatrixInterface.h"

namespace tudat
{

namespace propagators
{

//! Function to reset the state transition and sensitivity matrix interpolators
void SingleArcCombinedStateTransitionAndSensitivityMatrixInterface::updateMatrixInterpolators(
        const boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
        stateTransitionMatrixInterpolator,
        const boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
        sensitivityMatrixInterpolator )
{
    stateTransitionMatrixInterpolator_ = stateTransitionMatrixInterpolator;
    sensitivityMatrixInterpolator_ = sensitivityMatrixInterpolator;
}

//! Function to get the concatenated state transition and sensitivity matrix at a given time.
Eigen::MatrixXd SingleArcCombinedStateTransitionAndSensitivityMatrixInterface::getCombinedStateTransitionAndSensitivityMatrix(
        const double evaluationTime )
{
    combinedStateTransitionMatrix_.setZero( );

    // Set Phi and S matrices.
    combinedStateTransitionMatrix_.block( 0, 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_ ) =
            stateTransitionMatrixInterpolator_->interpolate( evaluationTime );

    if( sensitivityMatrixSize_ > 0 )
    {
        combinedStateTransitionMatrix_.block( 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_, sensitivityMatrixSize_ )=
                sensitivityMatrixInterpolator_->interpolate( evaluationTime );
    }

    return combinedStateTransitionMatrix_;
}

//! Constructor
MultiArcCombinedStateTransitionAndSensitivityMatrixInterface::MultiArcCombinedStateTransitionAndSensitivityMatrixInterface(
        const std::vector< boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
        stateTransitionMatrixInterpolators,
        const std::vector< boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
        sensitivityMatrixInterpolators,
        const std::vector< double >& arcStartTimes,
        const int numberOfInitialDynamicalParameters,
        const int numberOfParameters ):
    CombinedStateTransitionAndSensitivityMatrixInterface( numberOfInitialDynamicalParameters, numberOfParameters ),
    stateTransitionMatrixInterpolators_( stateTransitionMatrixInterpolators ),
    sensitivityMatrixInterpolators_( sensitivityMatrixInterpolators ),
    arcStartTimes_( arcStartTimes )
{
    numberOfStateArcs_ = arcStartTimes_.size( );

    sensitivityMatrixSize_ = numberOfParameters - numberOfStateArcs_ * stateTransitionMatrixSize_;

    if( stateTransitionMatrixInterpolators_.size( ) != sensitivityMatrixInterpolators_.size( ) ||
            stateTransitionMatrixInterpolators_.size( ) != static_cast< unsigned int >( numberOfStateArcs_ ) )
    {
        throw std::runtime_error(
                    "Error when making multi arc state transition and sensitivity interface, vector sizes are inconsistent" );
    }

    std::vector< double > arcSplitTimes = arcStartTimes_;
    arcSplitTimes.push_back(  std::numeric_limits< double >::max( ));
    lookUpscheme_ = boost::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >(
                arcSplitTimes );
}

//! Function to reset the state transition and sensitivity matrix interpolators
void MultiArcCombinedStateTransitionAndSensitivityMatrixInterface::updateMatrixInterpolators(
        const std::vector< boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
        stateTransitionMatrixInterpolators,
        const std::vector< boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
        sensitivityMatrixInterpolators,
        const std::vector< double >& arcStartTimes )
{
    stateTransitionMatrixInterpolators_ = stateTransitionMatrixInterpolators;
    sensitivityMatrixInterpolators_ = sensitivityMatrixInterpolators;
    arcStartTimes_ =  arcStartTimes;

    if( stateTransitionMatrixInterpolators_.size( ) != sensitivityMatrixInterpolators_.size( ) ||
            stateTransitionMatrixInterpolators_.size( ) != static_cast< unsigned int >( numberOfStateArcs_ ) )
    {
        throw std::runtime_error(
                    "Error when resetting multi arc state transition and sensitivity interface, vector sizes are inconsistent." );
    }

    std::vector< double > arcSplitTimes = arcStartTimes_;
    arcSplitTimes.push_back( std::numeric_limits< double >::max( ) );

    lookUpscheme_ = boost::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >(
                arcSplitTimes );
}

//! Function to get the concatenated single-arc state transition and sensitivity matrix at a given time.
Eigen::MatrixXd MultiArcCombinedStateTransitionAndSensitivityMatrixInterface::getCombinedStateTransitionAndSensitivityMatrix(
        const double evaluationTime )
{
    Eigen::MatrixXd combinedStateTransitionMatrix = Eigen::MatrixXd::Zero(
                stateTransitionMatrixSize_, stateTransitionMatrixSize_ + sensitivityMatrixSize_ );

    int currentArc = lookUpscheme_->findNearestLowerNeighbour( evaluationTime );

    // Set Phi and S matrices.
    combinedStateTransitionMatrix.block( 0, 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_ ) =
            stateTransitionMatrixInterpolators_.at( currentArc )->interpolate( evaluationTime );
    combinedStateTransitionMatrix.block( 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_, sensitivityMatrixSize_ ) =
            sensitivityMatrixInterpolators_.at( currentArc )->interpolate( evaluationTime );
    return combinedStateTransitionMatrix;
}

//! Function to get the concatenated state transition matrices for each arc and sensitivity matrix at a given time.
Eigen::MatrixXd MultiArcCombinedStateTransitionAndSensitivityMatrixInterface::getFullCombinedStateTransitionAndSensitivityMatrix(
        const double evaluationTime )
{
    Eigen::MatrixXd combinedStateTransitionMatrix = Eigen::MatrixXd::Zero(
                stateTransitionMatrixSize_, numberOfStateArcs_ * stateTransitionMatrixSize_ + sensitivityMatrixSize_ );

    int currentArc = lookUpscheme_->findNearestLowerNeighbour( evaluationTime );

    // Set Phi and S matrices of current arc.
    combinedStateTransitionMatrix.block( 0, currentArc * stateTransitionMatrixSize_, stateTransitionMatrixSize_, stateTransitionMatrixSize_ ) =
            stateTransitionMatrixInterpolators_.at( currentArc )->interpolate( evaluationTime );
    combinedStateTransitionMatrix.block(
                0, numberOfStateArcs_ * stateTransitionMatrixSize_, stateTransitionMatrixSize_, sensitivityMatrixSize_ ) =
            sensitivityMatrixInterpolators_.at( currentArc )->interpolate( evaluationTime );

    return combinedStateTransitionMatrix;
}

}

}

