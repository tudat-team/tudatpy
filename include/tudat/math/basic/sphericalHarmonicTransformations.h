#ifndef SPHERICALHARMONICTRANSFORMATIONS_H
#define SPHERICALHARMONICTRANSFORMATIONS_H

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <map>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/function.hpp>
#include <boost/math/special_functions/binomial.hpp>

#include "Tudat/Mathematics/BasicMathematics/wignerDMatrices.h"

namespace tudat
{

namespace basic_mathematics
{

class SphericalHarmonicTransformationCache
{
public:
    SphericalHarmonicTransformationCache(
            const int maxiumDegree,
            const int maximumOrder ):
        maximumDegree_( maxiumDegree ), maximumOrder_( maximumOrder ), updatePartials_( false )
    {
        wignerDMatricesCache_ = boost::make_shared< WignerDMatricesCache >( maxiumDegree );
    }

    void updateFromCayleyKleinParameters(
            const std::complex< double > cayleyKleinA,
            const std::complex< double > cayleyKleinB );

    void updateFromQuaternion(
            const Eigen::Quaterniond& rotationQuaternion );

    void updateFrom313EulerAngles(
            const double firstZRotation, const double XRotation, const double secondZRotation )
    {
        updateFromQuaternion(
                    Eigen::Quaterniond( Eigen::AngleAxisd( -secondZRotation, Eigen::Vector3d::UnitZ( ) ) *
                                        Eigen::AngleAxisd( -XRotation, Eigen::Vector3d::UnitX( ) ) *
                                        Eigen::AngleAxisd( -firstZRotation, Eigen::Vector3d::UnitZ( ) ) ) );
    }

    void transformCoefficientsAtDegree(
            const Eigen::MatrixXd& originalCosineCoefficients,
            const Eigen::MatrixXd& originalSineCoefficients,
            Eigen::MatrixXd& currentCosineCoefficients,
            Eigen::MatrixXd& currentSineCoefficients,
            const bool areCoefficientsNormalized = 1 );

    void setUpdatePartials( )
    {
        updatePartials_ = 1;
    }


private:

    std::complex< double > getWignerFunctionValue( const int l, const int m, const int k )
    {
        return wignerDMatricesCache_->getWignerDCoefficient( l, m, k );
    }

    boost::shared_ptr< WignerDMatricesCache > wignerDMatricesCache_;

    void setConstantCoefficients( );


    int maximumDegree_;

    int maximumOrder_;



    bool updatePartials_;

};

}

}

#endif // SPHERICALHARMONICTRANSFORMATIONS_H
