#ifndef SPHERICALHARMONICTRANSFORMATIONS_H
#define SPHERICALHARMONICTRANSFORMATIONS_H

#include <Eigen/Core>

#include <map>

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

    void transformCoefficientsAtDegree(
            const Eigen::MatrixXd& originalCosineCoefficients,
            const Eigen::MatrixXd& originalSineCoefficients,
            Eigen::MatrixXd& currentCosineCoefficients,
            Eigen::MatrixXd& currentSineCoefficients,
            const bool areCoefficientsNormalized = 1 );

    void getPartialDerivativesOfTransformedCoefficientsWrtEulerAngles(
            const Eigen::MatrixXd& originalCosineCoefficients,
            const Eigen::MatrixXd& originalSineCoefficients,
            std::vector< Eigen::MatrixXd >& currentCosineCoefficients,
            std::vector< Eigen::MatrixXd >& currentSineCoefficients,
            const bool areCoefficientsNormalized = 1 );

    void setUpdatePartials( )
    {
        updatePartials_ = 1;
    }


private:

    boost::shared_ptr< WignerDMatricesCache > wignerDMatricesCache_;

    void setConstantCoefficients( );


    int maximumDegree_;

    int maximumOrder_;



    bool updatePartials_;

};

}

}

#endif // SPHERICALHARMONICTRANSFORMATIONS_H
