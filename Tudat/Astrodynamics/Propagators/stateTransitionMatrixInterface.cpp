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

}

}

