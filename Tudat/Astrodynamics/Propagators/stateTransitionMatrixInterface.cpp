#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/Propagators/stateTransitionMatrixInterface.h"

namespace tudat
{

namespace propagators
{


    void SingleArcCombinedStateTransitionAndSensitivityMatrixInterface::updateMatrixInterpolators(
            const boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > stateTransitionMatrixInterpolator,
            const boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > sensitivityMatrixInterpolator )
    {
        stateTransitionMatrixInterpolator_ = stateTransitionMatrixInterpolator;
        sensitivityMatrixInterpolator_ = sensitivityMatrixInterpolator;
    }


    Eigen::MatrixXd SingleArcCombinedStateTransitionAndSensitivityMatrixInterface::getCombinesStateTransitionAndSensitivityMatrix(
            const double evaluationTime )
    {
        Eigen::MatrixXd combinedStateTransitionMatrix = Eigen::MatrixXd::Zero(
                    stateTransitionMatrixSize_, stateTransitionMatrixSize_ + sensitivityMatrixSize_ );

        // Set Phi and S matrices.
        combinedStateTransitionMatrix.block( 0, 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_ ) =
                stateTransitionMatrixInterpolator_->interpolate( evaluationTime );

        if( sensitivityMatrixSize_ > 0 )
        {
            combinedStateTransitionMatrix.block( 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_, sensitivityMatrixSize_ ) =
                    sensitivityMatrixInterpolator_->interpolate( evaluationTime );
        }

        return combinedStateTransitionMatrix;
    }

}

}

