/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_STATETRANSITIONMATRIXINTERFACE_H
#define TUDAT_STATETRANSITIONMATRIXINTERFACE_H

#include <iostream>
#include <vector>

#include <memory>

#include <Eigen/Core>

#include "tudat/math/interpolators/oneDimensionalInterpolator.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameterSet.h"

namespace tudat
{

namespace propagators
{

//! Base class for interface object of interpolation of numerically propagated state transition and sensitivity matrices.
/*!
 *  Base class for interface object of interpolation of numerically propagated state transition and sensitivity matrices.
 *  Derived classes implement the case of single-arc/multi-arc/combined dynamics.
 */
class CombinedStateTransitionAndSensitivityMatrixInterface
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param numberOfInitialDynamicalParameters Size of the estimated initial state vector (and size of square
     * state transition matrix.
     * \param numberOfParameters Total number of estimated parameters (initial states and other parameters).
     */
    CombinedStateTransitionAndSensitivityMatrixInterface(
            const int numberOfInitialDynamicalParameters,
            const int numberOfParameters )
    {
        stateTransitionMatrixSize_ = numberOfInitialDynamicalParameters;
        sensitivityMatrixSize_ = numberOfParameters - stateTransitionMatrixSize_;
    }

    //! Destructor.
    virtual ~CombinedStateTransitionAndSensitivityMatrixInterface( ){ }

    //! Function to get the concatenated state transition and sensitivity matrix at a given time.
    /*!
     *  Function to get the concatenated state transition and sensitivity matrix at a given time.
     *  \param evaluationTime Time at which to evaluate matrix interpolators
     *  \return Concatenated state transition and sensitivity matrices.
     */
    virtual Eigen::MatrixXd getCombinedStateTransitionAndSensitivityMatrix( const double evaluationTime,
                                                                            const std::vector< std::string >& arcDefiningBodies = std::vector< std::string >( ) ) = 0;

    //! Function to get the concatenated state transition and sensitivity matrix at a given time, which includes
    //! zero values for parameters not active in current arc.
    /*!
     *  Function to get the concatenated state transition and sensitivity matrix at a given time, which includes
     *  zero values for parameters not active in current arc.
     *  \param evaluationTime Time at which to evaluate matrix interpolators
     *  \return Concatenated state transition and sensitivity matrices, including inactive parameters at
     *  evaluationTime.
     */
    virtual Eigen::MatrixXd getFullCombinedStateTransitionAndSensitivityMatrix( const double evaluationTime,
                                                                                const std::vector< std::string >& arcDefiningBodies = std::vector< std::string >( ) ) = 0;

    //! Function to get the size of state transition matrix
    /*!
     * Function to get the size of state transition matrix
     * \return Size of state transition matrix
     */
    int getStateTransitionMatrixSize( )
    {
        return stateTransitionMatrixSize_;
    }

    //! Function to get the number of columns of sensitivity matrix.
    /*!
     * Function to get the number of columns of sensitivity matrix.
     * \return Number of columns of sensitivity matrix.
     */
    int getSensitivityMatrixSize( )
    {
        return sensitivityMatrixSize_;
    }

    //! Function to get the size of the total parameter vector.
    /*!
     * Function to get the size of the total parameter vector (both global and all local parameters).
     * \return Size of the total parameter vector.
     */
    virtual int getFullParameterVectorSize( ) = 0;

protected:

    //! Size of state transition matrix
    int stateTransitionMatrixSize_;

    //! Number of columns of sensitivity matrix.
    int sensitivityMatrixSize_;

};

//! Interface object of interpolation of numerically propagated state transition and sensitivity matrices for single-arc
//! estimation.
class SingleArcCombinedStateTransitionAndSensitivityMatrixInterface : public CombinedStateTransitionAndSensitivityMatrixInterface
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param stateTransitionMatrixInterpolator Interpolator returning the state transition matrix as a function of time.
     * \param sensitivityMatrixInterpolator Interpolator returning the sensitivity matrix as a function of time.
     * \param numberOfInitialDynamicalParameters Size of the estimated initial state vector (and size of square
     * state transition matrix.
     * \param numberOfParameters Total number of estimated parameters (initial states and other parameters).
     * \param statePartialAdditionIndices Vector of pair providing indices of column blocks of variational equations to add to other column blocks
     */
    SingleArcCombinedStateTransitionAndSensitivityMatrixInterface(
            const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
            stateTransitionMatrixInterpolator,
            const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
            sensitivityMatrixInterpolator,
            const int numberOfInitialDynamicalParameters,
            const int numberOfParameters,
            const std::vector< std::pair< int, int > >& statePartialAdditionIndices ):
        CombinedStateTransitionAndSensitivityMatrixInterface( numberOfInitialDynamicalParameters, numberOfParameters ),
        stateTransitionMatrixInterpolator_( stateTransitionMatrixInterpolator ),
        sensitivityMatrixInterpolator_( sensitivityMatrixInterpolator ),
        statePartialAdditionIndices_( statePartialAdditionIndices )
    {
        combinedStateTransitionMatrix_ = Eigen::MatrixXd::Zero(
                        stateTransitionMatrixSize_, stateTransitionMatrixSize_ + sensitivityMatrixSize_ );
    }

    //! Destructor.
    ~SingleArcCombinedStateTransitionAndSensitivityMatrixInterface( ){ }

    //! Function to reset the state transition and sensitivity matrix interpolators
    /*!
     * Function to reset the state transition and sensitivity matrix interpolators
     * \param stateTransitionMatrixInterpolator New interpolator returning the state transition matrix as a function of time.
     * \param sensitivityMatrixInterpolator New interpolator returning the sensitivity matrix as a function of time.
     * \param statePartialAdditionIndices Vector of pair providing indices of column blocks of variational equations to add to other column blocks
     */
    void updateMatrixInterpolators(
            const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
            stateTransitionMatrixInterpolator,
            const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
            sensitivityMatrixInterpolator,
            const std::vector< std::pair< int, int > >& statePartialAdditionIndices );

    //! Function to get the interpolator returning the state transition matrix as a function of time.
    /*!
     * Function to get the interpolator returning the state transition matrix as a function of time.
     * \return Interpolator returning the state transition matrix as a function of time.
     */
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
    getStateTransitionMatrixInterpolator( )
    {
        return stateTransitionMatrixInterpolator_;
    }

    //! Function to get the interpolator returning the sensitivity matrix as a function of time.
    /*!
     * Function to get the interpolator returning the sensitivity matrix as a function of time.
     * \return Interpolator returning the sensitivity matrix as a function of time.
     */
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
    getSensitivityMatrixInterpolator( )
    {
        return sensitivityMatrixInterpolator_;
    }

    //! Function to get the concatenated state transition and sensitivity matrix at a given time.
    /*!
     *  Function to get the concatenated state transition and sensitivity matrix at a given time.
     *  \param evaluationTime Time at which to evaluate matrix interpolators
     *  \return Concatenated state transition and sensitivity matrices.
     */
    Eigen::MatrixXd getCombinedStateTransitionAndSensitivityMatrix( const double evaluationTime,
                                                                    const std::vector< std::string >& arcDefiningBodies = std::vector< std::string >( ) );

    //! Function to get the concatenated state transition and sensitivity matrix at a given time.
    /*!
     *  Function to get the concatenated state transition and sensitivity matrix at a given time
     *  (functionality equal to getCombinedStateTransitionAndSensitivityMatrix for single-arc case).
     *  \param evaluationTime Time at which to evaluate matrix interpolators
     *  \return Concatenated state transition and sensitivity matrices.
     */
    Eigen::MatrixXd getFullCombinedStateTransitionAndSensitivityMatrix( const double evaluationTime,
                                                                        const std::vector< std::string >& arcDefiningBodies = std::vector< std::string >( ) )
    {
        return getCombinedStateTransitionAndSensitivityMatrix( evaluationTime, arcDefiningBodies );
    }

    //! Function to get the size of the total parameter vector.
    /*!
     * Function to get the size of the total parameter vector. For single-arc, this is simply the combination of
     * the size of the state transition and sensitivity matrices.
     * \return Size of the total parameter vector.
     */
    int getFullParameterVectorSize( )
    {
        return sensitivityMatrixSize_ + stateTransitionMatrixSize_;
    }

private:

    //! Predefined matrix to use as return value when calling getCombinedStateTransitionAndSensitivityMatrix.
    Eigen::MatrixXd combinedStateTransitionMatrix_;

    //! Interpolator returning the state transition matrix as a function of time.
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
    stateTransitionMatrixInterpolator_;

    //! Interpolator returning the sensitivity matrix as a function of time.
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
    sensitivityMatrixInterpolator_;

    std::vector< std::pair< int, int > > statePartialAdditionIndices_;

};

//! Interface object of interpolation of numerically propagated state transition and sensitivity matrices for multi-arc
//! estimation.
template< typename StateScalarType = double >
class MultiArcCombinedStateTransitionAndSensitivityMatrixInterface: public CombinedStateTransitionAndSensitivityMatrixInterface
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param stateTransitionMatrixInterpolators interpolators returning the state transition matrix as a function of time, vector
     * entries represent matrix history for each arc.
     * \param sensitivityMatrixInterpolators interpolators returning the sensitivity matrix as a function of time, vector
     * entries represent matrix history for each arc.
     * \param arcStartTimes Times at which the multiple arcs start
     * \param arcEndTimes Times at which the multiple arcs end
     * \param parametersToEstimate parameters to be estimated
     * \param numberOfInitialDynamicalParameters Size of the estimated initial state vector (and size of square
     * sing-arc state transition matrix times number of arcs.)
     * \param numberOfParameters Total number of estimated parameters (initial states and other parameters).
     * \param statePartialAdditionIndices Vector of pair providing indices of column blocks of variational equations to add to other column blocks
     */
    MultiArcCombinedStateTransitionAndSensitivityMatrixInterface(
            const std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
            stateTransitionMatrixInterpolators,
            const std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
            sensitivityMatrixInterpolators,
            const std::vector< double >& propagationStartTimes,
            const std::vector< double >& arcStartTimes,
            const std::vector< double >& arcEndTimes,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
            const int numberOfInitialDynamicalParameters,
            const int numberOfParameters,
            const std::vector< std::vector< std::pair< int, int > > >& statePartialAdditionIndices ):
        CombinedStateTransitionAndSensitivityMatrixInterface( numberOfInitialDynamicalParameters, numberOfParameters ),
        stateTransitionMatrixInterpolators_( stateTransitionMatrixInterpolators ),
        sensitivityMatrixInterpolators_( sensitivityMatrixInterpolators ),
        propagationStartTimes_( propagationStartTimes ),
        arcStartTimes_( arcStartTimes ),
        arcEndTimes_( arcEndTimes ),
        statePartialAdditionIndices_( statePartialAdditionIndices )
    {
        if( arcStartTimes_.size( ) != arcEndTimes_.size( ) )
        {
            throw std::runtime_error( "Error when making MultiArcCombinedStateTransitionAndSensitivityMatrixInterface, incompatible time lists" );
        }
        numberOfStateArcs_ = arcStartTimes_.size( );

        estimatable_parameters::getParametersToEstimatePerArcTest( parametersToEstimate, arcWiseParametersToEstimate_, propagationStartTimes_, estimatedBodiesPerArc_, arcIndicesPerBody_ );
        processArcWiseParametersIndices( parametersToEstimate, arcStartTimes_ );
        getArcStartTimesPerBody( );

        sensitivityMatrixSize_ = fullSensitivityMatrixSize_; // numberOfParameters - numberOfStateArcs_ * stateTransitionMatrixSize_;
        stateTransitionMatrixSize_ = fullStateSize_;
    //    std::cout << "TEST sensitivityMatrixSize_: " << sensitivityMatrixSize_ << "\n\n";
    //    std::cout << "TEST stateTransitionMatrixSize_: " << stateTransitionMatrixSize_ << "\n\n";

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

    //! Destructor
    ~MultiArcCombinedStateTransitionAndSensitivityMatrixInterface( ){ }

    //! Function to reset the state transition and sensitivity matrix interpolators
    /*!
     * Function to reset the state transition and sensitivity matrix interpolators
     * \param stateTransitionMatrixInterpolators New vector of interpolators returning the state transition matrix as a
     * function of time.
     * \param sensitivityMatrixInterpolators New vector of interpolator returning the sensitivity matrix as a function of time.
     * \param arcStartTimes Times at which the multiple arcs start.
     * \param arcEndTimes Times at which the multiple arcs end
     * \param statePartialAdditionIndices Vector of pair providing indices of column blocks of variational equations to add to other column blocks
     */
    void updateMatrixInterpolators(
            const std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
            stateTransitionMatrixInterpolators,
            const std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
            sensitivityMatrixInterpolators,
            const std::vector< double >& arcStartTimes,
            const std::vector< double >& arcEndTimes,
            const std::vector< std::vector< std::pair< int, int > > >& statePartialAdditionIndices  )
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

        lookUpschemePerBody_.clear( );
        for ( auto itr : arcStartTimesPerBody_ )
        {
            std::vector< double > arcSplitTimes = itr.second.first;
            arcSplitTimes.push_back(  std::numeric_limits< double >::max( ));
            lookUpschemePerBody_[ itr.first ] = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >( arcSplitTimes );
        }
    }

    //! Function to get the vector of interpolators returning the state transition matrix as a function of time.
    /*!
     * Function to get the vector of interpolators returning the state transition matrix as a function of time.
     * \return Vector of interpolators returning the state transition matrix as a function of time.
     */
    std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
    getStateTransitionMatrixInterpolators( )
    {
        return stateTransitionMatrixInterpolators_;
    }

    //! Function to get the vector of interpolators returning the sensitivity matrix as a function of time.
    /*!
     * Function to get the vector of interpolators returning the sensitivity matrix as a function of time.
     * \return Vector of interpolators returning the sensitivity matrix as a function of time.
     */
    std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
    getSensitivityMatrixInterpolator( )
    {
        return sensitivityMatrixInterpolators_;
    }

    //! Function to get the size of the total parameter vector.
    /*!
     * Function to get the size of the total parameter vector.
     * \return Size of the total parameter vector.
     */
    int getFullParameterVectorSize( )
    {
        return fullStateTransitionMatrixSize_ + fullSensitivityMatrixSize_; // sensitivityMatrixSize_ + numberOfStateArcs_ * stateTransitionMatrixSize_;
    }

    Eigen::MatrixXd getCombinedStateTransitionAndSensitivityMatrix(
            const double evaluationTime, const std::vector< std::string >& arcDefiningBodies = std::vector< std::string >( ),
            const bool addCentralBodySensitivity = true )
    {
        int currentArc = getCurrentArc( evaluationTime ).first;

        std::vector< int > currentArcsDefinedByEachBody;
        for ( unsigned int i = 0 ; i < arcDefiningBodies.size( ) ; i++ )
        {
            std::pair< int, double > currentArcDefinedByBody = getCurrentArc(evaluationTime, arcDefiningBodies.at( i ) );
    //        std::cout << "current arc defined by body: " << currentArcDefinedByBody.first << "\n\n";
            currentArcsDefinedByEachBody.push_back( currentArcDefinedByBody.first );
        }
        for ( unsigned int i = 0 ; i < currentArcsDefinedByEachBody.size( ) ; i++ )
        {
            if ( ( currentArcsDefinedByEachBody[ i ] != currentArcsDefinedByEachBody[ 0 ] ) && ( currentArcsDefinedByEachBody[ i ] != -1 )
                 && ( currentArcsDefinedByEachBody[ 0 ] != -1 ) )
            {
                std::runtime_error( "Error when getting current arc, different definitions for bodies " + arcDefiningBodies.at( i ) + " & " + arcDefiningBodies.at( 0 ) + "." );
            }
            if ( currentArcsDefinedByEachBody[ i ] != -1 )
            {
                currentArc = currentArcsDefinedByEachBody[ i ];
            }
        }

        int stateTransitionMatrixSize = 0;
        int sensitivityMatrixSize = 0;

        if ( currentArc >= 0 )
        {
            stateTransitionMatrixSize = arcWiseStateTransitionMatrixSize_[ currentArc ];
            sensitivityMatrixSize = arcWiseSensitivityMatrixSize_[ currentArc ];
        }

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
    //            std::cout << "multi-arc: statePartialAdditionIndices: " << statePartialAdditionIndices_.at( currentArc ).at( i ).first << " & " <<
    //                      statePartialAdditionIndices_.at( currentArc ).at( i ).second << " indicesToAdd: " << indicesToAdd << "\n\n";
                combinedStateTransitionMatrix.block(
                            statePartialAdditionIndices_.at( currentArc ).at( i ).first, 0, 6, indicesToAdd ) +=
                        combinedStateTransitionMatrix.block(
                            statePartialAdditionIndices_.at( currentArc ).at( i ).second, 0, 6, indicesToAdd );
            }
        }
        return combinedStateTransitionMatrix;
    }

    //! Function to get the concatenated single-arc state transition and sensitivity matrix at a given time.
    /*!
     *  Function to get the concatenated single-arc state transition and sensitivity matrix at a given time, evaluates matrices
     *  at the arc in which evaluationTime is located.
     *  \param evaluationTime Time at which to evaluate matrix interpolators
     *  \return Concatenated state transition and sensitivity matrices.
     */
    Eigen::MatrixXd getCombinedStateTransitionAndSensitivityMatrix(
            const double evaluationTime, const std::vector< std::string >& arcDefiningBodies = std::vector< std::string >( ) )
    {
        return getCombinedStateTransitionAndSensitivityMatrix( evaluationTime, arcDefiningBodies, true );
    }

    Eigen::MatrixXd getFullCombinedStateTransitionAndSensitivityMatrix(
            const double evaluationTime, const std::vector< std::string >& arcDefiningBodies = std::vector< std::string >( ),
            const bool addCentralBodySensitivity = true )
    {
        Eigen::MatrixXd combinedStateTransitionMatrix = getCombinedStateTransitionAndSensitivityMatrix(
                    evaluationTime, arcDefiningBodies, addCentralBodySensitivity );
        Eigen::MatrixXd fullCombinedStateTransitionMatrix = Eigen::MatrixXd::Zero(
    //                 stateTransitionMatrixSize_, numberOfStateArcs_ * stateTransitionMatrixSize_ + sensitivityMatrixSize_ );
                fullStateSize_, fullStateTransitionMatrixSize_ + fullSensitivityMatrixSize_ );

        int currentArc = getCurrentArc( evaluationTime ).first;

        std::vector< int > currentArcsDefinedByEachBody;
    //    std::cout << "arcDefiningBodies.size( ): " << arcDefiningBodies.size( ) << "\n\n";
        for ( unsigned int i = 0 ; i < arcDefiningBodies.size( ) ; i++ )
        {
            std::pair< int, double > currentArcDefinedByBody = getCurrentArc(evaluationTime, arcDefiningBodies.at( i ) );
    //        std::cout << "current arc defined by body: " << arcDefiningBodies.at( i ) << " = " << currentArcDefinedByBody.first << "\n\n";
            currentArcsDefinedByEachBody.push_back( currentArcDefinedByBody.first );
        }
        for ( unsigned int i = 0 ; i < currentArcsDefinedByEachBody.size( ) ; i++ )
        {
            if ( ( currentArcsDefinedByEachBody[ i ] != currentArcsDefinedByEachBody[ 0 ] ) && ( currentArcsDefinedByEachBody[ i ] != -1 )
                                                                                               && ( currentArcsDefinedByEachBody[ 0 ] != -1 ) )
            {
                throw std::runtime_error( "Error when getting current arc, different definitions for bodies " + arcDefiningBodies.at( i ) + " & " + arcDefiningBodies.at( 0 ) + "." );
            }
            if ( currentArcsDefinedByEachBody[ i ] != -1 )
            {
                currentArc = currentArcsDefinedByEachBody[ i ];
            }
        }

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
    /*!
     *  Function to get the concatenated state transition matrices for each arc and sensitivity matrix at a given time. The
     *  state transition matrix will be non-zero for only a single arc, but all state transition matrices at current arc are
     *  concatenated.
     *  \param evaluationTime Time at which to evaluate matrix interpolators
     *  \return Concatenated state transition and sensitivity matrices.
     */
    Eigen::MatrixXd getFullCombinedStateTransitionAndSensitivityMatrix( const double evaluationTime,
                                                                        const std::vector< std::string >& arcDefiningBodies = std::vector< std::string >( ) )
    {
        return getFullCombinedStateTransitionAndSensitivityMatrix( evaluationTime, arcDefiningBodies, true );
    }

    //! Function to retrieve the current arc for a given time
    /*!
     * Function to retrieve the current arc for a given time
     * \param evaluationTime Time at which current arc is to be determined
     * \return Pair with current arc index and associated arc initial time.
     */
    std::pair< int, double >  getCurrentArc( const double evaluationTime )
    {

        int currentArc =  lookUpscheme_->findNearestLowerNeighbour( evaluationTime );
        if( evaluationTime <= arcEndTimes_.at( currentArc ) && evaluationTime >= arcStartTimes_.at( currentArc ) )
        {
            return std::make_pair( currentArc, propagationStartTimes_.at( currentArc ) );
        }
        else
        {
            return std::make_pair( -1, TUDAT_NAN );
        }
    }

    //! Function to retrieve the current arc for a given time and a given body
    /*!
     * Function to retrieve the current arc for a given time and a given body
     * \param evaluationTime Time at which current arc is to be determined
     * \param body Name of the body for which current arc is to be determined
     * \return Pair with current arc index and associated arc initial time.
     */
    std::pair< int, double >  getCurrentArc( const double evaluationTime, const std::string body )
    {
        int currentArc = lookUpscheme_->findNearestLowerNeighbour( evaluationTime );
        if ( lookUpschemePerBody_.count( body ) != 0 )
        {
            currentArc = lookUpschemePerBody_.at( body )->findNearestLowerNeighbour( evaluationTime );
            if( evaluationTime < arcEndTimesPerBody_.at( body ).first.at( currentArc ) && evaluationTime > arcStartTimesPerBody_.at( body ).first.at( currentArc ) )
            {
    //            std::cout << "arcStartTimesPerBody_.at( body ).at( currentArc ): " << arcStartTimesPerBody_.at( body ).first.at( currentArc ) << "\n\n";
    //            std::cout << "arcEndTimesPerBody_.at( body ).at( currentArc ): " << arcEndTimesPerBody_.at( body ).first.at( currentArc ) << "\n\n";
    //            std::cout << "arcStartTimes_: " << arcStartTimes_.at( 1 ) << "\n\n";
    //            std::vector< double >::iterator itr = std::find( arcStartTimes_.begin( ), arcStartTimes_.end( ), arcStartTimesPerBody_.at( body ).at( currentArc ) );
    //            if ( itr != arcStartTimes_.cend( ) )
    //            {
    //                currentArc = std::distance( arcStartTimes_.begin( ), itr );
                    return std::make_pair( propagationStartTimesPerBody_.at( body ).second.at( currentArc ), propagationStartTimesPerBody_.at( body ).first.at( currentArc ) );
    //            }
    //            else
    //            {
    //                return std::make_pair( -1, TUDAT_NAN );
    //            }
            }
            else
            {
                return std::make_pair( -1, TUDAT_NAN );
            }
        }
        else
        {
    //        if( evaluationTime < arcEndTimes_.at( currentArc ) && evaluationTime > arcStartTimes_.at( currentArc ) )
    //        {
    //            return std::make_pair( currentArc, arcStartTimes_.at( currentArc ) );
    //        }
    //        else
    //        {
                return std::make_pair( -1, TUDAT_NAN );
    //        }
        }
    }

    //! Function to retrieve the number of arcs in dynamics
    /*!
     * Function to retrieve the number of arcs in dynamics
     * \return Number of arcs in dynamics
     */
    int getNumberOfArcs( )
    {
        return arcStartTimes_.size( );
    }

    std::vector< std::pair< int, int > > getStatePartialAdditionIndices(
            const int arcIndex )
    {
        return statePartialAdditionIndices_.at( arcIndex );
    }

    int getArcWiseStateTransitionMatrixSize( const int currentArc ) const
    {
        return arcWiseStateTransitionMatrixSize_[ currentArc ];
    }

    int getArcWiseSensitivityMatrixSize( const int currentArc ) const
    {
        return arcWiseSensitivityMatrixSize_[ currentArc ];
    }

    int getFullStateSize( ) const
    {
        return fullStateSize_;
    }

    int getFullStateTransitionMatrixSize( ) const
    {
        return fullStateTransitionMatrixSize_;
    }

    int getFullSensitivityMatrixSize( ) const
    {
        return fullSensitivityMatrixSize_;
    }

    std::map< int, std::map< std::string, std::pair< std::pair< int, int >, std::pair< std::pair< int, int >, int > > > > getArcWiseAndFullSolutionInitialStateIndices( ) const
    {
        return arcWiseAndFullSolutionInitialStateIndices_;
    }


protected:

    void processArcWiseParametersIndices( const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
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


            fullStateTransitionMatrixSize_ += arcWiseStateTransitionMatrixSize_[ arc ];

            ArcWiseParameterList estimatedBodiesPerArc = estimatable_parameters::getListOfBodiesWithTranslationalMultiArcStateToEstimate(
                    arcWiseParametersToEstimate_[ arc ] );
            for ( typename ArcWiseParameterList::const_iterator parameterIterator = estimatedBodiesPerArc.begin( ); parameterIterator !=
                                                                                                                    estimatedBodiesPerArc.end( ); parameterIterator++ )
            {
    //            std::cout << "estimated body: " << parameterIterator->first << "\n\n";
                estimatable_parameters::EstimatebleParameterIdentifier parameterId = parameterIterator->second->getParameterName( );
                std::vector< std::pair< int, int > > parametersIndices = arcWiseParametersToEstimate_[ arc ]->getIndicesForParameterType( parameterId );
                for ( unsigned int j = 0 ; j < parametersIndices.size( ) ; j++ )
                {
                    int startIndexInitialStateArcWise = parametersIndices[ j ].first;
                    int startIndexInitialStateFullSolution = indexBodyFullState[ parameterIterator->first ];
                    int sizeInitialState = parametersIndices[ j ].second;
                    int startIndexFullMatrix = initialStatesIndicesForFullMatrix.at( parameterIterator->first )[ 0 ].first
                                               + sizeInitialState * arcIndicesPerBody_.at( arc ).at( parameterIterator->first );
    //                std::cout << "in current arc: " << startIndexInitialStateArcWise << "\n\n";
    //                std::cout << "in full state: " << startIndexInitialStateFullSolution << "\n\n";
    //                std::cout << "sizeInitialState: " << sizeInitialState << "\n\n";
    //                std::cout << "in full matrix: " << startIndexFullMatrix << "\n\n";
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
    //    std::cout << "full sensitivity matrix size: " << fullSensitivityMatrixSize_ << "\n\n";
    //    std::cout << "full state transition matrix size: " << fullStateTransitionMatrixSize_ << "\n\n";
    //    std::cout << "full state size: " << fullStateSize_ << "\n\n";
    }

    void getArcStartTimesPerBody( )
    {
        for ( auto itr : estimatedBodiesPerArc_ )
        {
            for ( unsigned int i = 0 ; i < itr.second.size( ) ; i++ )
            {
    //            std::cout << "body: " << itr.second.at( i ) << "\n\n";
    //            std::cout << "arc start time per body: " << arcStartTimes_.at( itr.first ) << "\n\n";
                if ( arcStartTimesPerBody_.count( itr.second.at( i ) ) == 0 )
                {
                    std::vector< double > arcStartTimesVector = { arcStartTimes_.at( itr.first ) };
                    std::vector< double > arcEndTimesVector = { arcEndTimes_.at( itr.first ) };
                    std::vector< double > propagationStartTimesVector = { propagationStartTimes_.at( itr.first ) };
                    std::vector< int > arcIndicesVector = { itr.first };
                    arcStartTimesPerBody_[ itr.second.at( i ) ] = std::make_pair( arcStartTimesVector, arcIndicesVector );
                    arcEndTimesPerBody_[ itr.second.at( i ) ] = std::make_pair( arcEndTimesVector, arcIndicesVector );
                    propagationStartTimesPerBody_[ itr.second.at( i ) ] = std::make_pair( propagationStartTimesVector, arcIndicesVector );
                }
                else {
                    arcStartTimesPerBody_[ itr.second.at( i ) ].first.push_back( arcStartTimes_.at( itr.first ) );
                    arcStartTimesPerBody_[ itr.second.at( i) ].second.push_back( itr.first );
                    arcEndTimesPerBody_[ itr.second.at( i ) ].first.push_back( arcEndTimes_.at( itr.first ) );
                    arcEndTimesPerBody_[ itr.second.at( i ) ].second.push_back( itr.first );
                    propagationStartTimesPerBody_[ itr.second.at( i ) ].first.push_back( propagationStartTimes_.at( itr.first ) );
                    propagationStartTimesPerBody_[ itr.second.at( i) ].second.push_back( itr.first );
                }
            }
        }

        for ( auto itr : arcStartTimesPerBody_ )
        {
            std::vector< double > arcSplitTimes = itr.second.first;
            arcSplitTimes.push_back(  std::numeric_limits< double >::max( ));
            lookUpschemePerBody_[ itr.first ] = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >( arcSplitTimes );
        }
    }

private:

    //! List of interpolators returning the state transition matrix as a function of time.
    std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
    stateTransitionMatrixInterpolators_;

    //! List of interpolators returning the sensitivity matrix as a function of time.
    std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
    sensitivityMatrixInterpolators_;

    //! Times at which the propagation starts for each arc.
    std::vector< double > propagationStartTimes_;

    //! Times at which the multiple arcs start
    std::vector< double > arcStartTimes_;

    std::vector< double > arcEndTimes_;

    //! Number of arcs.
    int numberOfStateArcs_;

    //! Look-up algorithm to determine the arc of a given time.
    std::shared_ptr< interpolators::HuntingAlgorithmLookupScheme< double > > lookUpscheme_;

    //! Vector of pair providing indices of column blocks of variational equations to add to other column blocks
    /*!
     * Vector of pair providing indices of column blocks of variational equations to add to other column blocks,
     * which is needed by the hierarchical estimation fo dynamics. The second pair entry is the start of the column
     * block (of size 3) to where it should be copied. The first entry denotes from where it should be copied.
     * \sa setTranslationalStatePartialFrameScalingFunctions
     */
    std::vector< std::vector< std::pair< int, int > > > statePartialAdditionIndices_;

    //! Map containing, for each arc, a vector with the names of the bodies whose initial states are to be estimated.
    std::map< int, std::vector< std::string > > estimatedBodiesPerArc_;

    //! Vector containing the times at which each arc starts, per body.
    std::map< std::string, std::pair< std::vector< double >, std::vector< int > > > arcStartTimesPerBody_;

    //! Vector containing the times at which each arc ends, per body.
    std::map< std::string, std::pair< std::vector< double >, std::vector< int > > > arcEndTimesPerBody_;

    //! Vector containing the times at which each propagation starts, per body.
    std::map< std::string, std::pair< std::vector< double >, std::vector< int > > > propagationStartTimesPerBody_;

    //! Look-up algorithm to determine the arc of a given time, per body
    std::map< std::string, std::shared_ptr< interpolators::HuntingAlgorithmLookupScheme< double > > > lookUpschemePerBody_;

    //! map containing for each arc, the names of the bodies to be propagated/estimated and their arc-wise index for the given body.
    std::map< int, std::map< std::string, int > > arcIndicesPerBody_;

    //! size of the state transition matrix for current arc.
    std::vector< int > arcWiseStateTransitionMatrixSize_;

    //! size of the sensitivity matrix for current arc.
    std::vector< int > arcWiseSensitivityMatrixSize_;

    //! Vector of parameters to estimate, per arc.
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > > arcWiseParametersToEstimate_;

    //! map with indices of initial parameters in multi-arc solution and in the full combined solution
    //! key is the arc index. second element is a map: for each body to be estimated (key) in the given arc, we provide a pair of pairs, as follows:
    //! first pair: starting index & size of the initial state parameter in the given arc
    //! second pair: pair of starting index in full state and full matrix for the full combined solution & size of the initial state parameter
    std::map< int, std::map< std::string, std::pair< std::pair< int, int >, std::pair< std::pair< int, int >, int > > > > arcWiseAndFullSolutionInitialStateIndices_;

    //! size of the full initial state.
    int fullStateSize_;

    //! size of the full state transition matrix.
    int fullStateTransitionMatrixSize_;

    //! size of the full sensitivity matrix.
    int fullSensitivityMatrixSize_;

};

//! Interface object of interpolation of numerically propagated state transition and sensitivity matrices for a hybrid of
//! single-and multi-arc estimation (single order is put first in concatenation)
/*!
 *  Interface object of interpolation of numerically propagated state transition and sensitivity matrices for a hybrid of
 *  single-and multi-arc estimation (single order is put first in concatenation). The single- anf multi-arc given as input must
 *  be consistent: the single arc bodies/states must also be included in the multi-arc model, in order to properly generate
 *  the coupling terms between single and multi-arc states (see HybridArcVariationalEquationsSolver).
 */
template< typename StateScalarType = double >
class HybridArcCombinedStateTransitionAndSensitivityMatrixInterface: public CombinedStateTransitionAndSensitivityMatrixInterface
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param singleArcInterface Object to retrieve state transition/sensitivity matrices for single arc component
     * \param multiArcInterface Object to retrieve state transition/sensitivity matrices for multi arc component
     */
    HybridArcCombinedStateTransitionAndSensitivityMatrixInterface(
            const std::shared_ptr< SingleArcCombinedStateTransitionAndSensitivityMatrixInterface > singleArcInterface,
            const std::shared_ptr< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface< StateScalarType > > multiArcInterface ):
        CombinedStateTransitionAndSensitivityMatrixInterface(
            multiArcInterface->getStateTransitionMatrixSize( ),
            multiArcInterface->getSensitivityMatrixSize( ) + multiArcInterface->getStateTransitionMatrixSize( ) ),
        singleArcInterface_( singleArcInterface ), multiArcInterface_( multiArcInterface )
    {
        // Check input consistency
        if( multiArcInterface->getSensitivityMatrixSize( ) != singleArcInterface->getSensitivityMatrixSize( ) )
        {
            throw std::runtime_error( "Error when making hybrid state transition/sensitivity interface, input is inconsistent" );
        }
//        std::cout << "multi-arc sensitivity matrix size: " << multiArcInterface->getSensitivityMatrixSize( ) << " & single-arc sensitivity matrix size: "
//        << singleArcInterface->getSensitivityMatrixSize( ) << "\n\n";

        singleArcStateSize_ = singleArcInterface_->getStateTransitionMatrixSize( );
        multiArcStateSize_ = multiArcInterface_->getStateTransitionMatrixSize( );
        originalMultiArcStateSize_ = multiArcStateSize_ - singleArcStateSize_;
//        std::cout << "single arc state size: " << singleArcStateSize_ << "\n\n";
//        std::cout << "multi arc state size: " << multiArcStateSize_ << "\n\n";
//        std::cout << "original multi arc state size: " << originalMultiArcStateSize_ << "\n\n";

        numberOfMultiArcs_ = multiArcInterface->getNumberOfArcs( );


    }

    //! Destructor
    ~HybridArcCombinedStateTransitionAndSensitivityMatrixInterface( ){ }

    //! Function to get the size of the total parameter vector.
    /*!
     * Function to get the size of the total parameter vector.
     * \return Size of the total parameter vector.
     */
    int getFullParameterVectorSize( )
    {
//        std::cout << "multi arc - full parameter vector size: " << multiArcInterface_->getFullParameterVectorSize( ) << "\n\n";
//        std::cout << "hybrid arc - full parameter vector size: " << multiArcInterface_->getFullParameterVectorSize( ) + singleArcStateSize_ << "\n\n";
        return multiArcInterface_->getFullParameterVectorSize( ) - singleArcStateSize_ * ( numberOfMultiArcs_ - 1 );
        //sensitivityMatrixSize_ + singleArcStateSize_ + numberOfMultiArcs_ * originalMultiArcStateSize_;
    }

    //! Function to get the concatenated state transition and sensitivity matrix at a given time.
    /*!
     *  Function to get the concatenated state transition and sensitivity matrix at a given time. Only the state transition
     *  matrix for the current are is included in the concatenation.
     *  \param evaluationTime Time at which to evaluate matrix interpolators
     *  \return Concatenated state transition and sensitivity matrices.
     */
    Eigen::MatrixXd getCombinedStateTransitionAndSensitivityMatrix( const double evaluationTime,
                                                                    const std::vector< std::string >& arcDefiningBodies = std::vector< std::string >( ) )
    {
    //    std::cout << "before identifying current arc" << "\n\n";
        std::pair< int, double > currentArc = multiArcInterface_->getCurrentArc( evaluationTime );

        std::vector< std::pair< int, double > > currentArcsDefinedByEachBody;
        for ( unsigned int i = 0 ; i < arcDefiningBodies.size( ) ; i++ )
        {
            std::pair< int, double > currentArcDefinedByBody = multiArcInterface_->getCurrentArc(evaluationTime, arcDefiningBodies.at( i ) );
            currentArcsDefinedByEachBody.push_back( currentArcDefinedByBody );
        }
        for ( unsigned int i = 0 ; i < currentArcsDefinedByEachBody.size( ) ; i++ )
        {
            if ( ( currentArcsDefinedByEachBody[ i ] != currentArcsDefinedByEachBody[ 0 ] ) && ( currentArcsDefinedByEachBody[ i ].first != -1 )
                                                                                               && ( currentArcsDefinedByEachBody[ 0 ].first != -1 ) )
            {
                std::runtime_error( "Error when getting current arc, different definitions for bodies " + arcDefiningBodies.at( i ) + " & " + arcDefiningBodies.at( 0 ) + "." );
            }
            if ( currentArcsDefinedByEachBody[ i ].first != -1 )
            {
                currentArc.first = currentArcsDefinedByEachBody[ i ].first;
                currentArc.second = currentArcsDefinedByEachBody[ i ].second;
            }
        }


        int stateTransitionMatrixSize = singleArcInterface_->getStateTransitionMatrixSize( );
        int sensitivityMatrixSize = singleArcInterface_->getSensitivityMatrixSize( );
        if ( currentArc.first >= 0 )
        {
            stateTransitionMatrixSize = multiArcInterface_->getArcWiseStateTransitionMatrixSize( currentArc.first );
            sensitivityMatrixSize = multiArcInterface_->getArcWiseSensitivityMatrixSize( currentArc.first );
        }
        int multiArcStateSize = stateTransitionMatrixSize;
        int originalMultiArcStateSize = stateTransitionMatrixSize - singleArcStateSize_;
    //    std::cout << "multi-arc arc-wise state transition matrix size: " << stateTransitionMatrixSize << "\n\n";
    //    std::cout << "multi-arc arc-wise sensitivity matrix size: " << sensitivityMatrixSize << "\n\n";
    //    std::cout << "multi-arc arc-wise state size: " << multiArcStateSize << "\n\n";
    //    std::cout << "original multi-arc arc-wise state size: " << originalMultiArcStateSize << "\n\n";

        Eigen::MatrixXd combinedStateTransitionMatrix = Eigen::MatrixXd::Zero(
                    stateTransitionMatrixSize, stateTransitionMatrixSize + sensitivityMatrixSize );

        // Get single-arc matrices
    //    std::cout << "test0" << "\n\n";
        Eigen::MatrixXd singleArcStateTransition = singleArcInterface_->getCombinedStateTransitionAndSensitivityMatrix(
                    evaluationTime, arcDefiningBodies );

        // Set single-arc block
        combinedStateTransitionMatrix.block( 0, 0, singleArcStateSize_, singleArcStateSize_ ) =
                singleArcStateTransition.block( 0, 0, singleArcStateSize_, singleArcStateSize_ );

        // Set single-arc sensitivity block
        combinedStateTransitionMatrix.block( 0, multiArcStateSize, singleArcStateSize_, sensitivityMatrixSize_ ) =
                singleArcStateTransition.block( 0, singleArcStateSize_, singleArcStateSize_, sensitivityMatrixSize_ );

        if( !( currentArc.first < 0 ) )
        {
            // Get multi-arc matrices
    //        std::cout << "test1" << "\n\n";
            Eigen::MatrixXd multiArcStateTransition = multiArcInterface_->getCombinedStateTransitionAndSensitivityMatrix(
                    evaluationTime, arcDefiningBodies, false );

            // Set multi-arc block
            combinedStateTransitionMatrix.block(
                        singleArcStateSize_, singleArcStateSize_, originalMultiArcStateSize, originalMultiArcStateSize ) =
                    multiArcStateTransition.block(
                        singleArcStateSize_, singleArcStateSize_, originalMultiArcStateSize, originalMultiArcStateSize );

            // Get single-arc matrices at current arc start
    //        std::cout << "test2: " << currentArc.second << "\n\n";
            Eigen::MatrixXd singleArcStateTransitionAtArcStart = singleArcInterface_->getCombinedStateTransitionAndSensitivityMatrix(
                        currentArc.second, arcDefiningBodies );


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

    //        std::cout << "test3" << "\n\n";

            for( unsigned int i = 0; i < statePartialAdditionIndices.size( ); i++ )
            {
    //            std::cout << "state partial addition indices: " << statePartialAdditionIndices.at( i ).first << " & " << statePartialAdditionIndices.at( i ).second << "\n\n";
                combinedStateTransitionMatrix.block(
                            statePartialAdditionIndices.at( i ).first, multiArcStateSize,
                            6, sensitivityMatrixSize_ ) +=
                        combinedStateTransitionMatrix.block(
                            statePartialAdditionIndices.at( i ).second, multiArcStateSize,
                            6, sensitivityMatrixSize_ );

            }
    //        std::cout << "test4" << "\n\n";
        }

        return combinedStateTransitionMatrix;

    }

    //! Function to get the full concatenated state transition and sensitivity matrix at a given time.
    /*!
     *  Function to get the full concatenated state transition and sensitivity matrix at a given time. The state transition
     *  matrix for each arc is included (which equals zero for each multi-arc initial state sensitivity outside of teh current
     *  arc)
     *  \param evaluationTime Time at which to evaluate matrix interpolators
     *  \return Full concatenated state transition and sensitivity matrices.
     */
    Eigen::MatrixXd getFullCombinedStateTransitionAndSensitivityMatrix( const double evaluationTime,
                                                                        const std::vector< std::string >& arcDefiningBodies = std::vector< std::string >( ) )
    {
        int fullStateTransitionMatrixSize = multiArcInterface_->getFullStateTransitionMatrixSize( ) - singleArcStateSize_ * ( numberOfMultiArcs_ - 1 );
        int fullSensitivityMatrixSize = multiArcInterface_->getFullSensitivityMatrixSize( );

        Eigen::MatrixXd combinedStateTransitionMatrix = getCombinedStateTransitionAndSensitivityMatrix( evaluationTime, arcDefiningBodies );
        Eigen::MatrixXd fullCombinedStateTransitionMatrix = Eigen::MatrixXd::Zero(
                multiArcInterface_->getFullStateSize( ),
                multiArcInterface_->getFullStateTransitionMatrixSize( ) + multiArcInterface_->getFullSensitivityMatrixSize( ) - singleArcStateSize_ * ( numberOfMultiArcs_ - 1 ) );
    //    std::cout << "size full combined state transition matrix: " << fullCombinedStateTransitionMatrix.rows( ) << " & " << fullCombinedStateTransitionMatrix.cols( ) << "\n\n";
    //    std::cout << "sensitivity matrix size: " << sensitivityMatrixSize_ << "\n\n";
    //    std::cout << "multiArcInterface_->getFullSensitivityMatrixSize( ): " << multiArcInterface_->getFullSensitivityMatrixSize( ) << "\n\n";
    //    std::cout << "singleArcStateSize_: " << singleArcStateSize_ << "\n\n";
        std::pair< int, double > currentArc = multiArcInterface_->getCurrentArc( evaluationTime );

        std::vector< std::pair< int, double > > currentArcsDefinedByEachBody;
        for ( unsigned int i = 0 ; i < arcDefiningBodies.size( ) ; i++ )
        {
            std::pair< int, double > currentArcDefinedByBody = multiArcInterface_->getCurrentArc(evaluationTime, arcDefiningBodies.at( i ) );
    //        std::cout << "current arc defined by body: " << currentArcDefinedByBody.first << "\n\n";
            currentArcsDefinedByEachBody.push_back( currentArcDefinedByBody );
        }
        for ( unsigned int i = 0 ; i < currentArcsDefinedByEachBody.size( ) ; i++ )
        {
            if ( ( currentArcsDefinedByEachBody[ i ] != currentArcsDefinedByEachBody[ 0 ] ) && ( currentArcsDefinedByEachBody[ i ].first != -1 )
                                                                                               && ( currentArcsDefinedByEachBody[ 0 ].first != -1 ) )
            {
                throw std::runtime_error( "Error when getting current arc, different definitions for bodies " + arcDefiningBodies.at( i ) + " & " + arcDefiningBodies.at( 0 ) + "." );
            }
            if ( currentArcsDefinedByEachBody[ i ].first != -1  )
            {
                currentArc.first = currentArcsDefinedByEachBody[ i ].first;
                currentArc.second = currentArcsDefinedByEachBody[ i ].second;
            }
        }

        int stateTransitionMatrixSize = singleArcInterface_->getStateTransitionMatrixSize( );
        if ( currentArc.first >= 0 )
        {
            stateTransitionMatrixSize = multiArcInterface_->getArcWiseStateTransitionMatrixSize( currentArc.first );
        }
        int multiArcStateSize = stateTransitionMatrixSize;

        // Set single-arc block
        fullCombinedStateTransitionMatrix.block(
                    0, 0, singleArcStateSize_, singleArcStateSize_ ) =
                combinedStateTransitionMatrix.block(
                    0, 0, singleArcStateSize_, singleArcStateSize_ );

        // Set single-arc sensitivity block
        fullCombinedStateTransitionMatrix.block(
                    0, multiArcInterface_->getFullStateTransitionMatrixSize( ) - singleArcStateSize_ * ( numberOfMultiArcs_ - 1 ), singleArcStateSize_, sensitivityMatrixSize_ ) =
                combinedStateTransitionMatrix.block(
                    0, multiArcStateSize, singleArcStateSize_, sensitivityMatrixSize_ );

        // Set Phi and S matrices of current arc.
        if( currentArc.first >= 0 )
        {
            std::map< std::string, std::pair< std::pair< int, int >, std::pair< std::pair< int, int >, int > > > arcWiseAndFullSolutionIndices =
                    multiArcInterface_->getArcWiseAndFullSolutionInitialStateIndices( ).at( currentArc.first );

            // Set multi-arc block
            for ( auto itr : arcWiseAndFullSolutionIndices )
            {
    //            std::cout << "body: " << itr.first << "\n\n";
                std::pair< int, int > indicesInArcWiseSolution = itr.second.first;
                std::pair<std::pair< int, int >, int > indicesInFullSolution = itr.second.second;
                int indexInArcState = indicesInArcWiseSolution.first;
                int indexInArcMatrix = indicesInArcWiseSolution.second; // - singleArcStateSize_ * numberOfMultiArcs_;
                int indexInFullState = indicesInFullSolution.first.first;
                int indexInFullMatrix = indicesInFullSolution.first.second;
    //            std::cout << "indexInFullMatrix: " << indexInFullMatrix << "\n\n";
                int sizeInFullSolution = indicesInFullSolution.second;
    //            if ( indexInFullMatrix >= singleArcStateSize_ * ( currentArc.first + 1 ) )
    //            {
                    if ( indexInFullMatrix <= singleArcStateSize_ * numberOfMultiArcs_ )
                    {
                        indexInFullMatrix = indexInFullState;
                    }
                    else
                    {
                        indexInFullMatrix -= singleArcStateSize_ * ( numberOfMultiArcs_ - 1 );
                    }


    //                std::cout << "arc-wise indices: " << indexInArcState << " & " << indexInArcMatrix << "\n\n";
    //                std::cout << "full solution indices: " << indexInFullState << " & " << indexInFullMatrix << "\n\n";

                    // Set multi-arc block (self)
                    fullCombinedStateTransitionMatrix.block( indexInFullState, indexInFullMatrix, sizeInFullSolution, sizeInFullSolution ) =
                            combinedStateTransitionMatrix.block( indicesInArcWiseSolution.first, indicesInArcWiseSolution.first,
                                                                 indicesInArcWiseSolution.second, indicesInArcWiseSolution.second );

                    // Set coupled block
                    fullCombinedStateTransitionMatrix.block( indexInFullState, 0, sizeInFullSolution, sizeInFullSolution ) =
                            combinedStateTransitionMatrix.block( indicesInArcWiseSolution.first, 0, indicesInArcWiseSolution.second, indicesInArcWiseSolution.second );

                    // Set multi-arc sensitivity block
    //                std::cout << "sensitivity block: " << indexInFullState << " - " << fullStateTransitionMatrixSize << " & " <<
    //                indicesInFullSolution.second << " - " << fullSensitivityMatrixSize << "\n\n";
    //                std::cout << "sensitivity block: " << indicesInArcWiseSolution.first << " - " << multiArcInterface_->getArcWiseStateTransitionMatrixSize( currentArc.first ) << " & " <<
    //                indicesInArcWiseSolution.second << " - " << fullSensitivityMatrixSize << "\n\n";
                    fullCombinedStateTransitionMatrix.block(
                            indexInFullState, fullStateTransitionMatrixSize, indicesInFullSolution.second, fullSensitivityMatrixSize ) =
                            combinedStateTransitionMatrix.block( indicesInArcWiseSolution.first, multiArcInterface_->getArcWiseStateTransitionMatrixSize( currentArc.first ),
                                                                 indicesInArcWiseSolution.second, fullSensitivityMatrixSize );


                    // Set multi-arc block (other bodies)
                    for ( auto itr2 : arcWiseAndFullSolutionIndices )
                    {
                        if ( itr2.first != itr.first )
                        {
    //                        std::cout << "other body: " << itr2.first << "\n\n";
                            std::pair< int, int > indicesInArcWiseSolutionOtherBody = itr2.second.first;
                            std::pair< std::pair< int, int>, int > indicesInFullSolutionOtherBody = itr2.second.second;
                            int indexInArcStateOtherBody = indicesInArcWiseSolutionOtherBody.first;
                            int indexInFullStateOtherBody = indicesInFullSolutionOtherBody.first.first;
                            int indexInFullMatrixOtherBody = indicesInFullSolutionOtherBody.first.second;
                            int sizeInFullSolutionOtherBody = indicesInFullSolutionOtherBody.second;

    //                        if ( indexInFullMatrixOtherBody >= singleArcStateSize_ * ( currentArc.first + 1 ) )
    //                        {
                                if ( indexInFullMatrixOtherBody <= singleArcStateSize_ * numberOfMultiArcs_ )
                                {
                                    indexInFullMatrixOtherBody = indexInFullStateOtherBody;
                                }
                                else
                                {
                                    indexInFullMatrixOtherBody -= /*indexInArcStateOtherBody*/ singleArcStateSize_ * ( numberOfMultiArcs_ - 1 );
                                }

    //                            std::cout << "index full matrix other body: " << indexInFullMatrixOtherBody << "\n\n";

                            fullCombinedStateTransitionMatrix.block(
                                    indexInFullState, indexInFullMatrixOtherBody,
                                    indicesInFullSolution.second, sizeInFullSolutionOtherBody ) =
                                    combinedStateTransitionMatrix.block( indicesInArcWiseSolution.first, indicesInArcWiseSolutionOtherBody.first,
                                                                         indicesInArcWiseSolution.second, indicesInArcWiseSolutionOtherBody.second );

    //                        }
                        }
                    }
    //            }
            }

    //        // Set multi-arc block
    //        fullCombinedStateTransitionMatrix.block(
    //                    singleArcStateSize_, singleArcStateSize_ + currentArc.first * originalMultiArcStateSize_,
    //                    originalMultiArcStateSize_, originalMultiArcStateSize_ ) =
    //                combinedStateTransitionMatrix.block(
    //                    singleArcStateSize_, singleArcStateSize_, originalMultiArcStateSize_, originalMultiArcStateSize_ );
    //
    //
    //        // Set coupled block
    //        fullCombinedStateTransitionMatrix.block(
    //                    singleArcStateSize_, 0, originalMultiArcStateSize_, singleArcStateSize_ ) =
    //                combinedStateTransitionMatrix.block(
    //                    singleArcStateSize_, 0, originalMultiArcStateSize_, singleArcStateSize_ );
    //
    //        // Set multi-arc sensitivity block
    //        fullCombinedStateTransitionMatrix.block(
    //                    singleArcStateSize_, singleArcStateSize_ + numberOfMultiArcs_ * originalMultiArcStateSize_,
    //                    originalMultiArcStateSize_, sensitivityMatrixSize_ ) =
    //                combinedStateTransitionMatrix.block(
    //                    singleArcStateSize_, multiArcStateSize_, originalMultiArcStateSize_, sensitivityMatrixSize_ );

        }

        return fullCombinedStateTransitionMatrix;
    }


private:

    //! Object to retrieve state transition/sensitivity matrices for single arc component
    std::shared_ptr< SingleArcCombinedStateTransitionAndSensitivityMatrixInterface > singleArcInterface_;

    //! Object to retrieve state transition/sensitivity matrices for multi arc component
    std::shared_ptr< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface< StateScalarType > > multiArcInterface_;

    //! Size of single-arc state
    int singleArcStateSize_;

    //! Full state size of single arc in multi-arc model.
    int multiArcStateSize_;

    //! Full state size of single arc in original multi-arc model (full multi-arc size minus single-arc size).
    int originalMultiArcStateSize_;

    //! Number of arcs in multi-arc model.
    int numberOfMultiArcs_;
};



} // namespace propagators

} // namespace tudat

#endif // TUDAT_STATETRANSITIONMATRIXINTERFACE_H
