/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include "Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h"

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
    virtual Eigen::MatrixXd getCombinedStateTransitionAndSensitivityMatrix( const double evaluationTime ) = 0;

    //! Function to get the concatenated state transition and sensitivity matrix at a given time, which includes
    //! zero values for parameters not active in current arc.
    /*!
     *  Function to get the concatenated state transition and sensitivity matrix at a given time, which includes
     *  zero values for parameters not active in current arc.
     *  \param evaluationTime Time at which to evaluate matrix interpolators
     *  \return Concatenated state transition and sensitivity matrices, including inactive parameters at
     *  evaluationTime.
     */
    virtual Eigen::MatrixXd getFullCombinedStateTransitionAndSensitivityMatrix( const double evaluationTime ) = 0;

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
     */
    SingleArcCombinedStateTransitionAndSensitivityMatrixInterface(
            const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
            stateTransitionMatrixInterpolator,
            const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
            sensitivityMatrixInterpolator,
            const int numberOfInitialDynamicalParameters,
            const int numberOfParameters ):
        CombinedStateTransitionAndSensitivityMatrixInterface( numberOfInitialDynamicalParameters, numberOfParameters ),
        stateTransitionMatrixInterpolator_( stateTransitionMatrixInterpolator ),
        sensitivityMatrixInterpolator_( sensitivityMatrixInterpolator )
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
     */
    void updateMatrixInterpolators(
            const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
            stateTransitionMatrixInterpolator,
            const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
            sensitivityMatrixInterpolator );

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
    Eigen::MatrixXd getCombinedStateTransitionAndSensitivityMatrix( const double evaluationTime );

    //! Function to get the concatenated state transition and sensitivity matrix at a given time.
    /*!
     *  Function to get the concatenated state transition and sensitivity matrix at a given time
     *  (functionality equal to getCombinedStateTransitionAndSensitivityMatrix for single-arc case).
     *  \param evaluationTime Time at which to evaluate matrix interpolators
     *  \return Concatenated state transition and sensitivity matrices.
     */
    Eigen::MatrixXd getFullCombinedStateTransitionAndSensitivityMatrix( const double evaluationTime )
    {
        return getCombinedStateTransitionAndSensitivityMatrix( evaluationTime );
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
};

//! Interface object of interpolation of numerically propagated state transition and sensitivity matrices for multi-arc
//! estimation.
class MultiArcCombinedStateTransitionAndSensitivityMatrixInterface: public CombinedStateTransitionAndSensitivityMatrixInterface
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param stateTransitionMatrixInterpolators Interpolators returning the state transition matrix as a function of time, vector
     * entries represent matrix history for each arc.
     * \param sensitivityMatrixInterpolators Interpolators returning the sensitivity matrix as a function of time, vector
     * entries represent matrix history for each arc.
     * \param arcStartTimes Times at which the multiple arcs start
     * \param numberOfInitialDynamicalParameters Size of the estimated initial state vector (and size of square
     * sing-arc state transition matrix times number of arcs.)
     * \param numberOfParameters Total number of estimated parameters (initial states and other parameters).
     */
    MultiArcCombinedStateTransitionAndSensitivityMatrixInterface(
            const std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
            stateTransitionMatrixInterpolators,
            const std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
            sensitivityMatrixInterpolators,
            const std::vector< double >& arcStartTimes,
            const int numberOfInitialDynamicalParameters,
            const int numberOfParameters );

    //! Destructor
    ~MultiArcCombinedStateTransitionAndSensitivityMatrixInterface( ){ }

    //! Function to reset the state transition and sensitivity matrix interpolators
    /*!
     * Function to reset the state transition and sensitivity matrix interpolators
     * \param stateTransitionMatrixInterpolators New vector of interpolators returning the state transition matrix as a
     * function of time.
     * \param sensitivityMatrixInterpolators New vector of interpolator returning the sensitivity matrix as a function of time.
     * \param arcStartTimes Times at which the multiple arcs start.
     */
    void updateMatrixInterpolators(
            const std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
            stateTransitionMatrixInterpolators,
            const std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
            sensitivityMatrixInterpolators,
            const std::vector< double >& arcStartTimes );

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
        return sensitivityMatrixSize_ + numberOfStateArcs_ * stateTransitionMatrixSize_;
    }

    //! Function to get the concatenated single-arc state transition and sensitivity matrix at a given time.
    /*!
     *  Function to get the concatenated single-arc state transition and sensitivity matrix at a given time, evaluates matrices
     *  at the arc in which evaluationTime is located.
     *  \param evaluationTime Time at which to evaluate matrix interpolators
     *  \return Concatenated state transition and sensitivity matrices.
     */
    Eigen::MatrixXd getCombinedStateTransitionAndSensitivityMatrix( const double evaluationTime );

    //! Function to get the concatenated state transition matrices for each arc and sensitivity matrix at a given time.
    /*!
     *  Function to get the concatenated state transition matrices for each arc and sensitivity matrix at a given time. The
     *  state transition matrix will be non-zero for only a single arc, but all state transition matrices at current arc are
     *  concatenated.
     *  \param evaluationTime Time at which to evaluate matrix interpolators
     *  \return Concatenated state transition and sensitivity matrices.
     */
    Eigen::MatrixXd getFullCombinedStateTransitionAndSensitivityMatrix( const double evaluationTime );

    //! Function to retrieve the current arc for a given time
    /*!
     * Function to retrieve the current arc for a given time
     * \param evaluationTime Time at which current arc is to be determined
     * \return Pair with current arc index and associated arc initial time.
     */
    std::pair< int, double >  getCurrentArc( const double evaluationTime );

    //! Function to retrieve the number of arcs in dynamics
    /*!
     * Function to retrieve the number of arcs in dynamics
     * \return Number of arcs in dynamics
     */
    int getNumberOfArcs( )
    {
        return arcStartTimes_.size( );
    }

private:

    //! List of interpolators returning the state transition matrix as a function of time.
    std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
    stateTransitionMatrixInterpolators_;

    //! List of interpolators returning the sensitivity matrix as a function of time.
    std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
    sensitivityMatrixInterpolators_;

    //! Times at which the multiple arcs start
    std::vector< double > arcStartTimes_;

    //! Number of arcs.
    int numberOfStateArcs_;

    //! Look-up algorithm to determine the arc of a given time.
    std::shared_ptr< interpolators::HuntingAlgorithmLookupScheme< double > > lookUpscheme_;

};

//! Interface object of interpolation of numerically propagated state transition and sensitivity matrices for a hybrid of
//! single-and multi-arc estimation (single order is put first in concatenation)
/*!
 *  Interface object of interpolation of numerically propagated state transition and sensitivity matrices for a hybrid of
 *  single-and multi-arc estimation (single order is put first in concatenation). The single- anf multi-arc given as input must
 *  be consistent: the single arc bodies/states must also be included in the multi-arc model, in order to properly generate
 *  the coupling terms between single and multi-arc states (see HybridArcVariationalEquationsSolver).
 */
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
            const std::shared_ptr< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface > multiArcInterface ):
        CombinedStateTransitionAndSensitivityMatrixInterface(
            multiArcInterface->getStateTransitionMatrixSize( ),
            multiArcInterface->getSensitivityMatrixSize( ) + multiArcInterface->getStateTransitionMatrixSize( )),
        singleArcInterface_( singleArcInterface ), multiArcInterface_( multiArcInterface )
    {
        // Check input consistency
        if( multiArcInterface->getSensitivityMatrixSize( ) != singleArcInterface->getSensitivityMatrixSize( ) )
        {
            throw std::runtime_error( "Error when making hybrid state transition/sensitivity interface, input is inconsistent" );
        }

        singleArcStateSize_ = singleArcInterface_->getStateTransitionMatrixSize( );
        multiArcStateSize_ = multiArcInterface_->getStateTransitionMatrixSize( );
        originalMultiArcStateSize_ = multiArcStateSize_ - singleArcStateSize_;

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
        return sensitivityMatrixSize_ + singleArcStateSize_ + numberOfMultiArcs_ * originalMultiArcStateSize_;
    }

    //! Function to get the concatenated state transition and sensitivity matrix at a given time.
    /*!
     *  Function to get the concatenated state transition and sensitivity matrix at a given time. Only the state transition
     *  matrix for the current are is included in the concatenation.
     *  \param evaluationTime Time at which to evaluate matrix interpolators
     *  \return Concatenated state transition and sensitivity matrices.
     */
    Eigen::MatrixXd getCombinedStateTransitionAndSensitivityMatrix( const double evaluationTime );

    //! Function to get the full concatenated state transition and sensitivity matrix at a given time.
    /*!
     *  Function to get the full concatenated state transition and sensitivity matrix at a given time. The state transition
     *  matrix for each arc is included (which equals zero for each multi-arc initial state sensitivity outside of teh current
     *  arc)
     *  \param evaluationTime Time at which to evaluate matrix interpolators
     *  \return Full concatenated state transition and sensitivity matrices.
     */
    Eigen::MatrixXd getFullCombinedStateTransitionAndSensitivityMatrix( const double evaluationTime );

private:

    //! Object to retrieve state transition/sensitivity matrices for single arc component
    std::shared_ptr< SingleArcCombinedStateTransitionAndSensitivityMatrixInterface > singleArcInterface_;

    //! Object to retrieve state transition/sensitivity matrices for multi arc component
    std::shared_ptr< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface > multiArcInterface_;

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
