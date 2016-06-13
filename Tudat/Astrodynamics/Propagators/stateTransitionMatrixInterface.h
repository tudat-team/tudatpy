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


class CombinedStateTransitionAndSensitivityMatrixInterface
{
public:
    CombinedStateTransitionAndSensitivityMatrixInterface(
            const int numberOfInitialDynamicalParameters,
            const int numberOfParameters )
    {
        stateTransitionMatrixSize_ = numberOfInitialDynamicalParameters;
        sensitivityMatrixSize_ = numberOfParameters - stateTransitionMatrixSize_;
    }

    virtual ~CombinedStateTransitionAndSensitivityMatrixInterface( ){ }

    virtual Eigen::MatrixXd getCombinesStateTransitionAndSensitivityMatrix( const double evaluationTime ) = 0;

    //! Includes zeros for parameters not in current arc.
    virtual Eigen::MatrixXd getFullCombinesStateTransitionAndSensitivityMatrix( const double evaluationTime )
    {
        return getCombinesStateTransitionAndSensitivityMatrix( evaluationTime );
    }

    int getStateTransitionMatrixSize( )
    {
        return stateTransitionMatrixSize_;
    }

    int getSensitivityMatrixSize( )
    {
        return sensitivityMatrixSize_;
    }

    virtual int getFullParameterVectorSize( ) = 0;

protected:

    int stateTransitionMatrixSize_;

    int sensitivityMatrixSize_;

};

class SingleArcCombinedStateTransitionAndSensitivityMatrixInterface: public CombinedStateTransitionAndSensitivityMatrixInterface
{
public:
    SingleArcCombinedStateTransitionAndSensitivityMatrixInterface(
            const boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > stateTransitionMatrixInterpolator,
            const boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > sensitivityMatrixInterpolator,
            const int numberOfInitialDynamicalParameters,
            const int numberOfParameters ):
        CombinedStateTransitionAndSensitivityMatrixInterface( numberOfInitialDynamicalParameters, numberOfParameters ),
        stateTransitionMatrixInterpolator_( stateTransitionMatrixInterpolator ),
        sensitivityMatrixInterpolator_( sensitivityMatrixInterpolator )
    { }

    ~SingleArcCombinedStateTransitionAndSensitivityMatrixInterface( ){ }

    void updateMatrixInterpolators(
            const boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > stateTransitionMatrixInterpolator,
            const boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > sensitivityMatrixInterpolator );

    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
    getStateTransitionMatrixInterpolator( )
    {
        return stateTransitionMatrixInterpolator_;
    }

    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
    getSensitivityMatrixInterpolator( )
    {
        return sensitivityMatrixInterpolator_;
    }

    Eigen::MatrixXd getCombinesStateTransitionAndSensitivityMatrix( const double evaluationTime );

    int getFullParameterVectorSize( )
    {
        return sensitivityMatrixSize_ + stateTransitionMatrixSize_;
    }


private:

    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > stateTransitionMatrixInterpolator_;

    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > sensitivityMatrixInterpolator_;
};

}

}
#endif // STATETRANSITIONMATRIXINTERFACE_H
