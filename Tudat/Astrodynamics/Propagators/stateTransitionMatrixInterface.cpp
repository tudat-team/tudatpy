/*    Copyright (c) 2010-2017, Delft University of Technology
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
        const boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > stateTransitionMatrixInterpolator,
        const boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > sensitivityMatrixInterpolator )
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
        combinedStateTransitionMatrix_.block( 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_, sensitivityMatrixSize_ ) =
                sensitivityMatrixInterpolator_->interpolate( evaluationTime );
    }

    return combinedStateTransitionMatrix_;
}

MultiArcCombinedStateTransitionAndSensitivityMatrixInterface::MultiArcCombinedStateTransitionAndSensitivityMatrixInterface(
        const std::vector< boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > > stateTransitionMatrixInterpolators,
        const std::vector< boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > > sensitivityMatrixInterpolators,
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

    if( stateTransitionMatrixInterpolators_.size( ) != sensitivityMatrixInterpolators_.size( ) || stateTransitionMatrixInterpolators_.size( ) !=
            static_cast< unsigned int >( numberOfStateArcs_ ) )
    {
        std::cerr<<"Error when making multi arc state transition and sensitivity interface, vector sizes are inconsistent"<<std::endl;
    }

    std::vector< double > arcSplitTimes = arcStartTimes_;
    arcSplitTimes.push_back(  std::numeric_limits< double >::max( ));
    lookUpscheme_ = boost::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >(
                arcSplitTimes );
}

void MultiArcCombinedStateTransitionAndSensitivityMatrixInterface::updateMatrixInterpolators(
        const std::vector< boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > > stateTransitionMatrixInterpolator,
        const std::vector< boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > > sensitivityMatrixInterpolator,
        const std::vector< double >& arcStartTimes )
{
    stateTransitionMatrixInterpolators_ = stateTransitionMatrixInterpolator;
    sensitivityMatrixInterpolators_ = sensitivityMatrixInterpolator;
    arcStartTimes_ =  arcStartTimes;

    if( stateTransitionMatrixInterpolators_.size( ) != sensitivityMatrixInterpolators_.size( ) ||
            stateTransitionMatrixInterpolators_.size( ) != static_cast< unsigned int >( numberOfStateArcs_ ) )
    {
        std::cerr<<"Error when resetting multi arc state transition and sensitivity interface, vector sizes are inconsistent "<<
                   stateTransitionMatrixInterpolators_.size( )<<" "<<sensitivityMatrixInterpolators_.size( )<<" "<<
                   stateTransitionMatrixInterpolators_.size( )<<" "<<static_cast< unsigned int >( numberOfStateArcs_ )<<std::endl;
    }

    std::vector< double > arcSplitTimes = arcStartTimes_;
    arcSplitTimes.push_back( std::numeric_limits< double >::max( ) );

    lookUpscheme_ = boost::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >(
                arcSplitTimes );
}

Eigen::MatrixXd MultiArcCombinedStateTransitionAndSensitivityMatrixInterface::getCombinedStateTransitionAndSensitivityMatrix( const double evaluationTime )
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

Eigen::MatrixXd MultiArcCombinedStateTransitionAndSensitivityMatrixInterface::getFullCombinedStateTransitionAndSensitivityMatrix(
        const double evaluationTime )
{
    Eigen::MatrixXd combinedStateTransitionMatrix = Eigen::MatrixXd::Zero(
                stateTransitionMatrixSize_, numberOfStateArcs_ * stateTransitionMatrixSize_ + sensitivityMatrixSize_ );

    int currentArc = lookUpscheme_->findNearestLowerNeighbour( evaluationTime );

    combinedStateTransitionMatrix.block( 0, currentArc * stateTransitionMatrixSize_, stateTransitionMatrixSize_, stateTransitionMatrixSize_ ) =
            stateTransitionMatrixInterpolators_.at( currentArc )->interpolate( evaluationTime );

    combinedStateTransitionMatrix.block(
                0, numberOfStateArcs_ * stateTransitionMatrixSize_, stateTransitionMatrixSize_, sensitivityMatrixSize_ ) =
            sensitivityMatrixInterpolators_.at( currentArc )->interpolate( evaluationTime );

    return combinedStateTransitionMatrix;
}

std::pair< int, double > MultiArcCombinedStateTransitionAndSensitivityMatrixInterface::getCurrentArc( const double evaluationTime )
{
    int currentArc =  lookUpscheme_->findNearestLowerNeighbour( evaluationTime );
    return std::make_pair( currentArc, arcStartTimes_.at( currentArc ) );
}

Eigen::MatrixXd HybridArcCombinedStateTransitionAndSensitivityMatrixInterface::getCombinedStateTransitionAndSensitivityMatrix(
        const double evaluationTime )
{
    Eigen::MatrixXd combinedStateTransitionMatrix = Eigen::MatrixXd::Zero(
                stateTransitionMatrixSize_, stateTransitionMatrixSize_ + sensitivityMatrixSize_ );
    Eigen::MatrixXd singleArcStateTransition = singleArcInterface_->getCombinedStateTransitionAndSensitivityMatrix( evaluationTime );

    std::pair< int, double >  currentArc = multiArcInterface_->getCurrentArc( evaluationTime );

    Eigen::MatrixXd singleArcStateTransitionAtArcStart = singleArcInterface_->getCombinedStateTransitionAndSensitivityMatrix(
                currentArc.second );
    Eigen::MatrixXd multiArcStateTransition = multiArcInterface_->getCombinedStateTransitionAndSensitivityMatrix( evaluationTime );

    combinedStateTransitionMatrix.block( 0, 0, singleArcStateSize_, singleArcStateSize_ ) =
            singleArcStateTransition.block( 0, 0, singleArcStateSize_, singleArcStateSize_ );
    combinedStateTransitionMatrix.block(
                singleArcStateSize_, singleArcStateSize_, multiArcStateSize_ - singleArcStateSize_,
                multiArcStateSize_ - singleArcStateSize_ ) =
            multiArcStateTransition.block(
                singleArcStateSize_, singleArcStateSize_, multiArcStateSize_ - singleArcStateSize_,
                multiArcStateSize_ - singleArcStateSize_ );
    combinedStateTransitionMatrix.block(
                singleArcStateSize_, 0, multiArcStateSize_ - singleArcStateSize_,
                singleArcStateSize_ ) =
            multiArcStateTransition.block(
                singleArcStateSize_, 0, multiArcStateSize_ - singleArcStateSize_,
                singleArcStateSize_ ) * singleArcStateTransitionAtArcStart.block(
                0, 0, singleArcStateSize_, singleArcStateSize_ );
    combinedStateTransitionMatrix.block( 0, multiArcStateSize_, singleArcStateSize_, sensitivityMatrixSize_ ) =
            singleArcStateTransition.block( 0, singleArcStateSize_, singleArcStateSize_, sensitivityMatrixSize_ );
    combinedStateTransitionMatrix.block(
                singleArcStateSize_, multiArcStateSize_, multiArcStateSize_ - singleArcStateSize_, sensitivityMatrixSize_ ) =
            multiArcStateTransition.block(
                singleArcStateSize_, multiArcStateSize_, multiArcStateSize_ - singleArcStateSize_, sensitivityMatrixSize_ );

    return combinedStateTransitionMatrix;

}

Eigen::MatrixXd HybridArcCombinedStateTransitionAndSensitivityMatrixInterface::getFullCombinedStateTransitionAndSensitivityMatrix(
        const double evaluationTime )
{
    Eigen::MatrixXd combinedStateTransitionMatrix = Eigen::MatrixXd::Zero(
                stateTransitionMatrixSize_, singleArcStateSize_ + numberOfMultiArcs_ * multiArcStateSize_ + sensitivityMatrixSize_ );

    std::pair< int, double >  currentArc = multiArcInterface_->getCurrentArc( evaluationTime );

    Eigen::MatrixXd singleArcStateTransitionSensitivity =
            singleArcInterface_->getCombinedStateTransitionAndSensitivityMatrix( evaluationTime );
    Eigen::MatrixXd singleArcStateTransitionSensitivityAtArcStart =
            singleArcInterface_->getCombinedStateTransitionAndSensitivityMatrix( currentArc.second );
    Eigen::MatrixXd multiArcStateTransitionSensitivity =
            multiArcInterface_->getCombinedStateTransitionAndSensitivityMatrix( evaluationTime );

    // Add single-arc bodies own dependencies
    combinedStateTransitionMatrix.block( 0, 0, singleArcStateSize_, singleArcStateSize_ ) =
            singleArcStateTransitionSensitivity;

    // Add multi-arc bodies own dependencies
    combinedStateTransitionMatrix.block(
                singleArcStateSize_, singleArcStateSize_ + currentArc.first * multiArcStateSize_, multiArcStateSize_, multiArcStateSize_ ) =
            multiArcStateTransitionSensitivity.block(
                singleArcStateSize_, singleArcStateSize_, multiArcStateSize_, multiArcStateSize_ );

    // Add multi-arc bodies own dependencies
    combinedStateTransitionMatrix.block(
                singleArcStateSize_, 0, multiArcStateSize_, singleArcStateSize_ ) =
            multiArcStateTransitionSensitivity.block( singleArcStateSize_, 0, multiArcStateSize_, singleArcStateSize_ ) *
            singleArcStateTransitionSensitivityAtArcStart;

    // Add sensitivity matrix components
    combinedStateTransitionMatrix.block(
                0, singleArcStateSize_ + numberOfMultiArcs_ * multiArcStateSize_, singleArcStateSize_, sensitivityMatrixSize_ ) =
          singleArcStateTransitionSensitivity.block( 0, singleArcStateSize_, singleArcStateSize_, sensitivityMatrixSize_ );
    combinedStateTransitionMatrix.block(
                singleArcStateSize_, singleArcStateSize_ + numberOfMultiArcs_ * multiArcStateSize_, multiArcStateSize_, sensitivityMatrixSize_ ) =
          singleArcStateTransitionSensitivity.block(
                singleArcStateSize_, multiArcStateSize_, multiArcStateSize_, sensitivityMatrixSize_ );

    return combinedStateTransitionMatrix;
}

}

}

