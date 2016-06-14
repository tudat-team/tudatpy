#ifndef STATETRANSITIONMATRIXINTERFACE_H
#define STATETRANSITIONMATRIXINTERFACE_H

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
     * \brief Function to get the interpolator returning the state transition matrix as a function of time.
     * \return Interpolator returning the state transition matrix as a function of time.
     */
    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
    getStateTransitionMatrixInterpolator( )
    {
        return stateTransitionMatrixInterpolator_;
    }

    //! Function to get the interpolator returning the sensitivity matrix as a function of time.
    /*!
     * \brief Function to get the interpolator returning the sensitivity matrix as a function of time.
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

}

}
#endif // STATETRANSITIONMATRIXINTERFACE_H
