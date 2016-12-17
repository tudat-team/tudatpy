#ifndef PODINPUTOUTPUTTYPES_H
#define PODINPUTOUTPUTTYPES_H

#include <map>
#include <vector>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/ObservationModels/observableTypes.h"

namespace tudat
{

namespace simulation_setup
{

template< typename ObservationScalarType = double, typename TimeType = double, typename StateScalarType = ObservationScalarType,
          typename ParameterScalarType = double >
struct PodInput
{

public:
    typedef Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > ObservationVectorType;
    typedef std::map< observation_models::LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< TimeType >, observation_models::LinkEndType > > > SingleObservablePodInputType;
    typedef std::map< observation_models::ObservableType, SingleObservablePodInputType > PodInputDataType;

    PodInput( const PodInputDataType& observationsAndTimes,
              const Eigen::Matrix< ParameterScalarType, Eigen::Dynamic, 1 > initialParameterDeviationEstimate,
              const bool reintegrateEquationsOnFirstIteration = 1,
              const Eigen::MatrixXd inverseOfAprioriCovariance = Eigen::MatrixXd::Zero( 0, 0 ),
              const double constantWeight = 1.0 ):
        observationsAndTimes_( observationsAndTimes ), initialParameterDeviationEstimate_( initialParameterDeviationEstimate ),
        reintegrateEquationsOnFirstIteration_( reintegrateEquationsOnFirstIteration ), inverseOfAprioriCovariance_( inverseOfAprioriCovariance )
    {
        if( inverseOfAprioriCovariance_.rows( ) == 0 )
        {
            setDefaultAprioriInverseCovariance( );
        }
        setDefaultWeightsMatrix( constantWeight );
    }


    void setDefaultAprioriInverseCovariance(  )
    {
        inverseOfAprioriCovariance_ = Eigen::MatrixXd::Constant(
                    initialParameterDeviationEstimate_.rows( ), initialParameterDeviationEstimate_.rows( ), 0.0 );
    }

    void setDefaultWeightsMatrix( const double constantWeight = 0.0 )
    {
        for( typename PodInputDataType::const_iterator observablesIterator = observationsAndTimes_.begin( );
             observablesIterator != observationsAndTimes_.end( ); observablesIterator++ )
        {
            for( typename SingleObservablePodInputType::const_iterator dataIterator =
                 observablesIterator->second.begin( ); dataIterator != observablesIterator->second.end( ); dataIterator++  )
            {
                weightsMatrixDiagonals_[ observablesIterator->first ][ dataIterator->first ] =
                        Eigen::VectorXd::Constant( dataIterator->second.first.rows( ), constantWeight );
            }
        }
    }

    PodInputDataType observationsAndTimes_;

    Eigen::Matrix< ParameterScalarType, Eigen::Dynamic, 1 > initialParameterDeviationEstimate_;

    bool reintegrateEquationsOnFirstIteration_;

    Eigen::MatrixXd inverseOfAprioriCovariance_;

    std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, Eigen::VectorXd > > weightsMatrixDiagonals_;


};


template< typename ParameterScalarType = double >
struct PodOutput
{
    PodOutput( const Eigen::Matrix< ParameterScalarType, Eigen::Dynamic, 1 >& parameterEstimate,
               const Eigen::VectorXd& residuals,
               const Eigen::MatrixXd& normalizedInformationMatrix,
               const Eigen::VectorXd& weightsMatrixDiagonal,
               const Eigen::VectorXd& informationMatrixTransformationDiagonal,
               const Eigen::MatrixXd& inverseNormalizedCovarianceMatrix,
               const double meanAbsoluteResidual ):

        parameterEstimate_( parameterEstimate ), residuals_( residuals ),
        normalizedInformationMatrix_( normalizedInformationMatrix ), weightsMatrixDiagonal_( weightsMatrixDiagonal ),
        informationMatrixTransformationDiagonal_( informationMatrixTransformationDiagonal ),
        inverseNormalizedCovarianceMatrix_( inverseNormalizedCovarianceMatrix ),
        meanAbsoluteResidual_( meanAbsoluteResidual ){ }

    Eigen::Matrix< ParameterScalarType, Eigen::Dynamic, 1 > parameterEstimate_;
    Eigen::VectorXd residuals_;
    Eigen::MatrixXd normalizedInformationMatrix_;
    Eigen::VectorXd weightsMatrixDiagonal_;
    Eigen::VectorXd informationMatrixTransformationDiagonal_;
    Eigen::MatrixXd inverseNormalizedCovarianceMatrix_;
    double meanAbsoluteResidual_;
};

}

}
#endif // PODINPUTOUTPUTTYPES_H
