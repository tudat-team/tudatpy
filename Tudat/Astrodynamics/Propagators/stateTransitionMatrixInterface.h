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

#include <vector>

#include <boost/shared_ptr.hpp>

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
class SingleArcCombinedStateTransitionAndSensitivityMatrixInterface:
        public CombinedStateTransitionAndSensitivityMatrixInterface
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
            const boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
            stateTransitionMatrixInterpolator,
            const boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
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
            const boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
            stateTransitionMatrixInterpolator,
            const boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
            sensitivityMatrixInterpolator );

    //! Function to get the interpolator returning the state transition matrix as a function of time.
    /*!
     * Function to get the interpolator returning the state transition matrix as a function of time.
     * \return Interpolator returning the state transition matrix as a function of time.
     */
    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
    getStateTransitionMatrixInterpolator( )
    {
        return stateTransitionMatrixInterpolator_;
    }

    //! Function to get the interpolator returning the sensitivity matrix as a function of time.
    /*!
     * Function to get the interpolator returning the sensitivity matrix as a function of time.
     * \return Interpolator returning the sensitivity matrix as a function of time.
     */
    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
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
    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
    stateTransitionMatrixInterpolator_;

    //! Interpolator returning the sensitivity matrix as a function of time.
    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
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
            const std::vector< boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
            stateTransitionMatrixInterpolators,
            const std::vector< boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
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
            const std::vector< boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
            stateTransitionMatrixInterpolators,
            const std::vector< boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
            sensitivityMatrixInterpolators,
            const std::vector< double >& arcStartTimes );

    //! Function to get the vector of interpolators returning the state transition matrix as a function of time.
    /*!
     * Function to get the vector of interpolators returning the state transition matrix as a function of time.
     * \return Vector of interpolators returning the state transition matrix as a function of time.
     */
    std::vector< boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
    getStateTransitionMatrixInterpolators( )
    {
        return stateTransitionMatrixInterpolators_;
    }

    //! Function to get the vector of interpolators returning the sensitivity matrix as a function of time.
    /*!
     * Function to get the vector of interpolators returning the sensitivity matrix as a function of time.
     * \return Vector of interpolators returning the sensitivity matrix as a function of time.
     */
    std::vector< boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
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

private:

    //! List of interpolators returning the state transition matrix as a function of time.
    std::vector< boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
    stateTransitionMatrixInterpolators_;

    //! List of interpolators returning the sensitivity matrix as a function of time.
    std::vector< boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
    sensitivityMatrixInterpolators_;

    //! Times at which the multiple arcs start
    std::vector< double > arcStartTimes_;

    //! Number of arcs.
    int numberOfStateArcs_;

    //! Look-up algorithm to determine the arc of a given time.
    boost::shared_ptr< interpolators::HuntingAlgorithmLookupScheme< double > > lookUpscheme_;

};


} // namespace propagators

} // namespace tudat
#endif // TUDAT_STATETRANSITIONMATRIXINTERFACE_H
