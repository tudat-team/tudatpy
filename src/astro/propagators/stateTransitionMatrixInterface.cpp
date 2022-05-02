/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/make_shared.hpp>

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
    statePartialAdditionIndices_ = statePartialAdditionIndices;
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

    for( unsigned int i = 0; i < statePartialAdditionIndices_.size( ); i++ )
    {
        combinedStateTransitionMatrix_.block(
                    statePartialAdditionIndices_.at( i ).first, 0, 6, stateTransitionMatrixSize_ + sensitivityMatrixSize_ ) +=
                combinedStateTransitionMatrix_.block(
                    statePartialAdditionIndices_.at( i ).second, 0, 6, stateTransitionMatrixSize_ + sensitivityMatrixSize_ );
    }


    return combinedStateTransitionMatrix_;
}

//! Constructor
MultiArcCombinedStateTransitionAndSensitivityMatrixInterface::MultiArcCombinedStateTransitionAndSensitivityMatrixInterface(
        const std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
        stateTransitionMatrixInterpolators,
        const std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
        sensitivityMatrixInterpolators,
        const std::vector< double >& arcStartTimes,
        const std::vector< double >& arcEndTimes,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate,
        const int numberOfInitialDynamicalParameters,
        const int numberOfParameters,
        const std::vector< std::vector< std::pair< int, int > > >& statePartialAdditionIndices ):
    CombinedStateTransitionAndSensitivityMatrixInterface( numberOfInitialDynamicalParameters, numberOfParameters ),
    stateTransitionMatrixInterpolators_( stateTransitionMatrixInterpolators ),
    sensitivityMatrixInterpolators_( sensitivityMatrixInterpolators ),
    arcStartTimes_( arcStartTimes ),
    arcEndTimes_( arcEndTimes ),
    statePartialAdditionIndices_( statePartialAdditionIndices )
{
    if( arcStartTimes_.size( ) != arcEndTimes_.size( ) )
    {
        throw std::runtime_error( "Error when making MultiArcCombinedStateTransitionAndSensitivityMatrixInterface, incompatible time lists" );
    }
    numberOfStateArcs_ = arcStartTimes_.size( );

    getParametersToEstimatePerArcTest( parametersToEstimate, arcWiseParametersToEstimate_, arcStartTimes_, estimatedBodiesPerArc_, arcIndicesPerBody_ );
    processArcWiseParametersIndices( parametersToEstimate, arcStartTimes_ );

    sensitivityMatrixSize_ = fullSensitivityMatrixSize_; // numberOfParameters - numberOfStateArcs_ * stateTransitionMatrixSize_;
    stateTransitionMatrixSize_ = fullStateSize_;
    std::cout << "TEST sensitivityMatrixSize_: " << sensitivityMatrixSize_ << "\n\n";
    std::cout << "TEST stateTransitionMatrixSize_: " << stateTransitionMatrixSize_ << "\n\n";

    if( stateTransitionMatrixInterpolators_.size( ) != sensitivityMatrixInterpolators_.size( ) ||
            stateTransitionMatrixInterpolators_.size( ) != static_cast< unsigned int >( numberOfStateArcs_ ) )
    {
        throw std::runtime_error(
                    "Error when making multi arc state transition and sensitivity interface, vector sizes are inconsistent" );
    }

    std::vector< double > arcSplitTimes = arcStartTimes_;
    arcSplitTimes.push_back(  std::numeric_limits< double >::max( ));
    lookUpscheme_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >(
                arcSplitTimes );
}

//! Function to reset the state transition and sensitivity matrix interpolators
void MultiArcCombinedStateTransitionAndSensitivityMatrixInterface::updateMatrixInterpolators(
        const std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
        stateTransitionMatrixInterpolators,
        const std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
        sensitivityMatrixInterpolators,
        const std::vector< double >& arcStartTimes,
        const std::vector< double >& arcEndTimes,
        const std::vector< std::vector< std::pair< int, int > > >& statePartialAdditionIndices )
{
    stateTransitionMatrixInterpolators_ = stateTransitionMatrixInterpolators;
    sensitivityMatrixInterpolators_ = sensitivityMatrixInterpolators;
    arcStartTimes_ =  arcStartTimes;
    arcEndTimes_ = arcEndTimes;
    statePartialAdditionIndices_ = statePartialAdditionIndices;

    if( stateTransitionMatrixInterpolators_.size( ) != sensitivityMatrixInterpolators_.size( ) ||
            stateTransitionMatrixInterpolators_.size( ) != static_cast< unsigned int >( numberOfStateArcs_ ) )
    {
        throw std::runtime_error(
                    "Error when resetting multi arc state transition and sensitivity interface, vector sizes are inconsistent." );
    }

    std::vector< double > arcSplitTimes = arcStartTimes_;
    arcSplitTimes.push_back( std::numeric_limits< double >::max( ) );

    lookUpscheme_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >(
                arcSplitTimes );
}

//! Function to get the concatenated single-arc state transition and sensitivity matrix at a given time.
Eigen::MatrixXd MultiArcCombinedStateTransitionAndSensitivityMatrixInterface::getCombinedStateTransitionAndSensitivityMatrix(
        const double evaluationTime,
        const bool addCentralBodySensitivity )
{
    int currentArc = getCurrentArc( evaluationTime ).first;

    int stateTransitionMatrixSize = arcWiseStateTransitionMatrixSize_[ currentArc ];
    int sensitivityMatrixSize = arcWiseSensitivityMatrixSize_[ currentArc ];

    Eigen::MatrixXd combinedStateTransitionMatrix = Eigen::MatrixXd::Zero(
                stateTransitionMatrixSize, stateTransitionMatrixSize + sensitivityMatrixSize );



    // Set Phi and S matrices.
    if( currentArc >= 0 )
    {
        combinedStateTransitionMatrix.block( 0, 0, stateTransitionMatrixSize, stateTransitionMatrixSize ) =
                stateTransitionMatrixInterpolators_.at( currentArc )->interpolate( evaluationTime );
        combinedStateTransitionMatrix.block( 0, stateTransitionMatrixSize, stateTransitionMatrixSize, sensitivityMatrixSize ) =
                sensitivityMatrixInterpolators_.at( currentArc )->interpolate( evaluationTime );

        for( unsigned int i = 0; i < statePartialAdditionIndices_.at( currentArc ).size( ); i++ )
        {
            int indicesToAdd = addCentralBodySensitivity ?
                        ( stateTransitionMatrixSize + sensitivityMatrixSize ) : stateTransitionMatrixSize;
            combinedStateTransitionMatrix.block(
                        statePartialAdditionIndices_.at( currentArc ).at( i ).first, 0, 6, indicesToAdd ) +=
                    combinedStateTransitionMatrix.block(
                        statePartialAdditionIndices_.at( currentArc ).at( i ).second, 0, 6, indicesToAdd );
        }
    }
    return combinedStateTransitionMatrix;
}

//! Function to get the concatenated single-arc state transition and sensitivity matrix at a given time.
Eigen::MatrixXd MultiArcCombinedStateTransitionAndSensitivityMatrixInterface::getCombinedStateTransitionAndSensitivityMatrix(
        const double evaluationTime )
{
    return getCombinedStateTransitionAndSensitivityMatrix( evaluationTime, true );
}

//! Function to get the concatenated state transition matrices for each arc and sensitivity matrix at a given time.
Eigen::MatrixXd MultiArcCombinedStateTransitionAndSensitivityMatrixInterface::getFullCombinedStateTransitionAndSensitivityMatrix(
        const double evaluationTime,
        const bool addCentralBodySensitivity )
{
    Eigen::MatrixXd combinedStateTransitionMatrix = getCombinedStateTransitionAndSensitivityMatrix(
                evaluationTime, addCentralBodySensitivity );
    Eigen::MatrixXd fullCombinedStateTransitionMatrix = Eigen::MatrixXd::Zero(
//                 stateTransitionMatrixSize_, numberOfStateArcs_ * stateTransitionMatrixSize_ + sensitivityMatrixSize_ );
            fullStateSize_, fullStateTransitionMatrixSize_ + fullSensitivityMatrixSize_ );

    int currentArc = getCurrentArc( evaluationTime ).first;

    // Set Phi and S matrices of current arc.
    if( currentArc >= 0 )
    {
//        fullCombinedStateTransitionMatrix.block(
//                    0, currentArc * stateTransitionMatrixSize_, stateTransitionMatrixSize_, stateTransitionMatrixSize_ ) =
//                combinedStateTransitionMatrix.block( 0, 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_ );
//        fullCombinedStateTransitionMatrix.block(
//                    0, numberOfStateArcs_ * stateTransitionMatrixSize_, stateTransitionMatrixSize_, sensitivityMatrixSize_ ) =
//                combinedStateTransitionMatrix.block( 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_, sensitivityMatrixSize_ );

        std::map< std::string, std::pair< std::pair< int, int >, std::pair< std::pair< int, int >, int > > > arcWiseAndFullSolutionIndices =
                arcWiseAndFullSolutionInitialStateIndices_.at( currentArc );
        for ( auto itr : arcWiseAndFullSolutionIndices )
        {
            std::pair< int, int > indicesInArcWiseSolution = itr.second.first;
            std::pair< std::pair< int, int >, int > indicesInFullSolution = itr.second.second;
            int indexInFullState = indicesInFullSolution.first.first;
            int indexInFullMatrix = indicesInFullSolution.first.second;
            int sizeInFullSolution = indicesInFullSolution.second;

            fullCombinedStateTransitionMatrix.block(
                    indexInFullState, indexInFullMatrix,
                    sizeInFullSolution, sizeInFullSolution ) =
                    combinedStateTransitionMatrix.block( indicesInArcWiseSolution.first, indicesInArcWiseSolution.first,
                                                         indicesInArcWiseSolution.second, indicesInArcWiseSolution.second );

            for ( auto itr2 : arcWiseAndFullSolutionIndices )
            {
                if ( itr2.first != itr.first )
                {
                    std::pair< int, int > indicesInArcWiseSolutionOtherBody = itr2.second.first;
                    std::pair< std::pair< int, int>, int > indicesInFullSolutionOtherBody = itr2.second.second;
//                    int indexInFullStateOtherBody = indicesInFullSolutionOtherBody.first.first;
                    int indexInFullMatrixOtherBody = indicesInFullSolutionOtherBody.first.second;
                    int sizeInFullSolutionOtherBody = indicesInFullSolutionOtherBody.second;

                    fullCombinedStateTransitionMatrix.block(
                            indexInFullState, indexInFullMatrixOtherBody,
                            indicesInFullSolution.second, sizeInFullSolutionOtherBody ) =
                            combinedStateTransitionMatrix.block( indicesInArcWiseSolution.first, indicesInArcWiseSolutionOtherBody.first,
                                                                 indicesInArcWiseSolution.second, indicesInArcWiseSolutionOtherBody.second );
                }
            }

            fullCombinedStateTransitionMatrix.block(
                    indexInFullState, fullStateTransitionMatrixSize_, indicesInFullSolution.second, fullSensitivityMatrixSize_ ) =
                    combinedStateTransitionMatrix.block( indicesInArcWiseSolution.first, arcWiseStateTransitionMatrixSize_[ currentArc ],
                                                         indicesInArcWiseSolution.second, arcWiseSensitivityMatrixSize_[ currentArc ] );
        }


    }
    return fullCombinedStateTransitionMatrix;
}

//! Function to get the concatenated state transition matrices for each arc and sensitivity matrix at a given time.
Eigen::MatrixXd MultiArcCombinedStateTransitionAndSensitivityMatrixInterface::getFullCombinedStateTransitionAndSensitivityMatrix(
        const double evaluationTime )
{
    return getFullCombinedStateTransitionAndSensitivityMatrix( evaluationTime, true );
}

//! Function to retrieve the current arc for a given time
std::pair< int, double > MultiArcCombinedStateTransitionAndSensitivityMatrixInterface::getCurrentArc( const double evaluationTime )
{

    int currentArc =  lookUpscheme_->findNearestLowerNeighbour( evaluationTime );
    if( evaluationTime < arcEndTimes_.at( currentArc ) && evaluationTime > arcStartTimes_.at( currentArc ) )
    {
        return std::make_pair( currentArc, arcStartTimes_.at( currentArc ) );
    }
    else
    {
        return std::make_pair( -1, TUDAT_NAN );
    }
}


template< typename StateScalarType, typename TimeType >
void MultiArcCombinedStateTransitionAndSensitivityMatrixInterface::processArcWiseParametersIndices(
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
        const std::vector< double > arcStartTimes )
{

    // Get list of objets and associated bodies to estimate initial arc-wise translational states
    typedef std::map< std::string, std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > > ArcWiseParameterList;
    ArcWiseParameterList estimatedBodies = estimatable_parameters::getListOfBodiesWithTranslationalMultiArcStateToEstimate(
            parametersToEstimate );

    // Iterate over all parameters
    unsigned int counterEstimatedBody = 0;
    fullStateSize_ = 0;
    std::map< std::string, int > indexBodyFullState;
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > > multiArcInitialStateParameters
    = parametersToEstimate->getEstimatedMultiArcInitialStateParameters( );
    for ( unsigned int i = 0 ; i < multiArcInitialStateParameters.size( ) ; i++ )
    {
        indexBodyFullState[ multiArcInitialStateParameters[ i ]->getParameterName( ).second.first ] = counterEstimatedBody * 6;
        counterEstimatedBody ++;
        fullStateSize_ += 6;
    }

    std::map< std::string, std::vector< std::pair< int, int > > > initialStatesIndicesForFullMatrix;
    for ( typename ArcWiseParameterList::const_iterator parameterIterator = estimatedBodies.begin( );
    parameterIterator != estimatedBodies.end( ); parameterIterator++ )
    {
//        std::cout << "estimated body: " << parameterIterator->first << "\n\n";
        estimatable_parameters::EstimatebleParameterIdentifier parameterId = parameterIterator->second->getParameterName( );
        std::vector< std::pair< int, int > > parametersIndices = parametersToEstimate->getIndicesForParameterType( parameterId );
        if (  parametersIndices.size( ) != 1 )
        {
            throw std::runtime_error( "Error when creating state transition matrix interface, more than one parameter found with ID " +
            parameterIterator->second->getParameterDescription( ) );
        }

        initialStatesIndicesForFullMatrix[ parameterIterator->first ] = parametersIndices;
    }

    fullStateTransitionMatrixSize_ = 0;
    fullSensitivityMatrixSize_ = 0;
    for ( unsigned int arc = 0 ; arc < estimatedBodiesPerArc_.size( ) ; arc++ )
    {
        arcWiseStateTransitionMatrixSize_.push_back( getSingleArcInitialDynamicalStateParameterSetSize( arcWiseParametersToEstimate_[ arc ] ) );
        arcWiseSensitivityMatrixSize_.push_back( getSingleArcParameterSetSize( arcWiseParametersToEstimate_[ arc ] ) - arcWiseStateTransitionMatrixSize_[ arc ] );

        std::cout << "arc " << arc << " - arcWiseStateTransitionMatrixSize_: " << arcWiseStateTransitionMatrixSize_[ arc ] << "\n\n";
        std::cout << "arc " << arc << " - arcWiseSensitivityMatrixSize_: " << arcWiseSensitivityMatrixSize_[ arc ] << "\n\n";

        fullStateTransitionMatrixSize_ += arcWiseStateTransitionMatrixSize_[ arc ];

        ArcWiseParameterList estimatedBodiesPerArc = estimatable_parameters::getListOfBodiesWithTranslationalMultiArcStateToEstimate(
                arcWiseParametersToEstimate_[ arc ] );
        for ( typename ArcWiseParameterList::const_iterator parameterIterator = estimatedBodiesPerArc.begin( ); parameterIterator !=
                                                                                                                estimatedBodiesPerArc.end( ); parameterIterator++ )
        {
            std::cout << "estimated body: " << parameterIterator->first << "\n\n";
            estimatable_parameters::EstimatebleParameterIdentifier parameterId = parameterIterator->second->getParameterName( );
            std::vector< std::pair< int, int > > parametersIndices = arcWiseParametersToEstimate_[ arc ]->getIndicesForParameterType( parameterId );
            for ( unsigned int j = 0 ; j < parametersIndices.size( ) ; j++ )
            {
                int startIndexInitialStateArcWise = parametersIndices[ j ].first;
                int startIndexInitialStateFullSolution = indexBodyFullState[ parameterIterator->first ];
                int sizeInitialState = parametersIndices[ j ].second;
                int startIndexFullMatrix = initialStatesIndicesForFullMatrix.at( parameterIterator->first )[ 0 ].first
                                           + sizeInitialState * arcIndicesPerBody_.at( arc ).at( parameterIterator->first );
                std::cout << "in current arc: " << startIndexInitialStateArcWise << "\n\n";
                std::cout << "in full state: " << startIndexInitialStateFullSolution << "\n\n";
                std::cout << "sizeInitialState: " << sizeInitialState << "\n\n";
                std::cout << "in full matrix: " << startIndexFullMatrix << "\n\n";
                arcWiseAndFullSolutionInitialStateIndices_[ arc ][ parameterIterator->first ] =
                        std::make_pair( std::make_pair( startIndexInitialStateArcWise, sizeInitialState ),
                                        std::make_pair( std::make_pair( startIndexInitialStateFullSolution, startIndexFullMatrix), sizeInitialState ) );
            }
        }
    }

    for ( unsigned int arc = 1 ; arc < estimatedBodiesPerArc_.size( ) ; arc++ )
    {
        if ( arcWiseSensitivityMatrixSize_[ arc ] != arcWiseSensitivityMatrixSize_[ 0 ] )
        {
            throw std::runtime_error( "Current implementation cannot (yet) handle varying sensitivity matrix size from one arc to another.");
        }
    }
    fullSensitivityMatrixSize_ = arcWiseSensitivityMatrixSize_[ 0 ];
    std::cout << "full sensitivity matrix size: " << fullSensitivityMatrixSize_ << "\n\n";
    std::cout << "full state transition matrix size: " << fullStateTransitionMatrixSize_ << "\n\n";
    std::cout << "full state size: " << fullStateSize_ << "\n\n";
}


//! Function to get the concatenated state transition and sensitivity matrix at a given time.
Eigen::MatrixXd HybridArcCombinedStateTransitionAndSensitivityMatrixInterface::getCombinedStateTransitionAndSensitivityMatrix(
        const double evaluationTime )
{
    std::pair< int, double > currentArc = multiArcInterface_->getCurrentArc( evaluationTime );
    int stateTransitionMatrixSize = multiArcInterface_->getArcWiseStateTransitionMatrixSize( currentArc.first );
    int sensitivityMatrixSize = multiArcInterface_->getArcWiseSensitivityMatrixSize( currentArc.first );
    int multiArcStateSize = stateTransitionMatrixSize;
    int originalMultiArcStateSize = stateTransitionMatrixSize - singleArcStateSize_;
    std::cout << "multi-arc arc-wise state transition matrix size: " << stateTransitionMatrixSize << "\n\n";
    std::cout << "multi-arc arc-wise sensitivity matrix size: " << sensitivityMatrixSize << "\n\n";
    std::cout << "multi-arc arc-wise state size: " << multiArcStateSize << "\n\n";
    std::cout << "original multi-arc arc-wise state size: " << originalMultiArcStateSize << "\n\n";

    Eigen::MatrixXd combinedStateTransitionMatrix = Eigen::MatrixXd::Zero(
                stateTransitionMatrixSize, stateTransitionMatrixSize + sensitivityMatrixSize );

    // Get single-arc matrices
    Eigen::MatrixXd singleArcStateTransition = singleArcInterface_->getCombinedStateTransitionAndSensitivityMatrix(
                evaluationTime );

    // Get multi-arc matrices
    Eigen::MatrixXd multiArcStateTransition = multiArcInterface_->getCombinedStateTransitionAndSensitivityMatrix(
                evaluationTime, false );

    // Set single-arc block
    combinedStateTransitionMatrix.block( 0, 0, singleArcStateSize_, singleArcStateSize_ ) =
            singleArcStateTransition.block( 0, 0, singleArcStateSize_, singleArcStateSize_ );

    // Set single-arc sensitivity block
    combinedStateTransitionMatrix.block( 0, multiArcStateSize, singleArcStateSize_, sensitivityMatrixSize_ ) =
            singleArcStateTransition.block( 0, singleArcStateSize_, singleArcStateSize_, sensitivityMatrixSize_ );

    if( !( currentArc.first < 0 ) )
    {
        // Set multi-arc block
        combinedStateTransitionMatrix.block(
                    singleArcStateSize_, singleArcStateSize_, originalMultiArcStateSize, originalMultiArcStateSize ) =
                multiArcStateTransition.block(
                    singleArcStateSize_, singleArcStateSize_, originalMultiArcStateSize, originalMultiArcStateSize );

        // Get single-arc matrices at current arc start
        Eigen::MatrixXd singleArcStateTransitionAtArcStart = singleArcInterface_->getCombinedStateTransitionAndSensitivityMatrix(
                    currentArc.second );

        // Set coupled block
        combinedStateTransitionMatrix.block(
                    singleArcStateSize_, 0, originalMultiArcStateSize, singleArcStateSize_ ) =
                multiArcStateTransition.block(
                    singleArcStateSize_, 0, originalMultiArcStateSize, singleArcStateSize_ ) *
                singleArcStateTransitionAtArcStart.block(
                    0, 0, singleArcStateSize_, singleArcStateSize_ );

        // Set multi-arc sensitivity block
        combinedStateTransitionMatrix.block(
                    singleArcStateSize_, multiArcStateSize, originalMultiArcStateSize, sensitivityMatrixSize_ ) =
                multiArcStateTransition.block(
                    singleArcStateSize_, multiArcStateSize, originalMultiArcStateSize, sensitivityMatrixSize_ );

        std::vector< std::pair< int, int > > statePartialAdditionIndices =
                multiArcInterface_->getStatePartialAdditionIndices( currentArc.first );

        for( unsigned int i = 0; i < statePartialAdditionIndices.size( ); i++ )
        {
            std::cout << "state partial addition indices: " << statePartialAdditionIndices.at( i ).first << " & " << statePartialAdditionIndices.at( i ).second << "\n\n";
            combinedStateTransitionMatrix.block(
                        statePartialAdditionIndices.at( i ).first, multiArcStateSize,
                        6, sensitivityMatrixSize_ ) +=
                    combinedStateTransitionMatrix.block(
                        statePartialAdditionIndices.at( i ).second, multiArcStateSize,
                        6, sensitivityMatrixSize_ );

        }
    }

    return combinedStateTransitionMatrix;

}

//! Function to get the concatenated state transition matrices for each arc and sensitivity matrix at a given time.
Eigen::MatrixXd HybridArcCombinedStateTransitionAndSensitivityMatrixInterface::getFullCombinedStateTransitionAndSensitivityMatrix(
        const double evaluationTime )
{
    Eigen::MatrixXd combinedStateTransitionMatrix = getCombinedStateTransitionAndSensitivityMatrix( evaluationTime );
    Eigen::MatrixXd fullCombinedStateTransitionMatrix = Eigen::MatrixXd::Zero(
                stateTransitionMatrixSize_, singleArcStateSize_ + numberOfMultiArcs_ * originalMultiArcStateSize_ + sensitivityMatrixSize_ );

    fullCombinedStateTransitionMatrix.block(
                0, 0, singleArcStateSize_, singleArcStateSize_ ) =
            combinedStateTransitionMatrix.block(
                0, 0, singleArcStateSize_, singleArcStateSize_ );

    fullCombinedStateTransitionMatrix.block(
                0, singleArcStateSize_ + numberOfMultiArcs_ * originalMultiArcStateSize_, singleArcStateSize_, sensitivityMatrixSize_ ) =
            combinedStateTransitionMatrix.block(
                0, multiArcStateSize_, singleArcStateSize_, sensitivityMatrixSize_ );

    int currentArc = multiArcInterface_->getCurrentArc( evaluationTime ).first;

    // Set Phi and S matrices of current arc.
    if( currentArc >= 0 )
    {

        // Set multi-arc block
        fullCombinedStateTransitionMatrix.block(
                    singleArcStateSize_, singleArcStateSize_ + currentArc * originalMultiArcStateSize_,
                    originalMultiArcStateSize_, originalMultiArcStateSize_ ) =
                combinedStateTransitionMatrix.block(
                    singleArcStateSize_, singleArcStateSize_, originalMultiArcStateSize_, originalMultiArcStateSize_ );


        // Set coupled block
        fullCombinedStateTransitionMatrix.block(
                    singleArcStateSize_, 0, originalMultiArcStateSize_, singleArcStateSize_ ) =
                combinedStateTransitionMatrix.block(
                    singleArcStateSize_, 0, originalMultiArcStateSize_, singleArcStateSize_ );

        // Set multi-arc sensitivity block
        fullCombinedStateTransitionMatrix.block(
                    singleArcStateSize_, singleArcStateSize_ + numberOfMultiArcs_ * originalMultiArcStateSize_,
                    originalMultiArcStateSize_, sensitivityMatrixSize_ ) =
                combinedStateTransitionMatrix.block(
                    singleArcStateSize_, multiArcStateSize_, originalMultiArcStateSize_, sensitivityMatrixSize_ );

    }

    return fullCombinedStateTransitionMatrix;
}

}

}

