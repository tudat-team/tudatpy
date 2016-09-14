#ifndef SPHERICALHARMONICPARTIALFUNCTIONS_H
#define SPHERICALHARMONICPARTIALFUNCTIONS_H

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/sphericalHarmonics.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{

void computePotentialSphericalHessian(
        const double distance,
        const double radiusPowerTerm,
        const double cosineOfOrderLongitude,
        const double sineOfOrderLongitude,
        const double cosineOfLatitude,
        const double sineOfLatitude,
        const double preMultiplier,
        const int degree,
        const int order,
        const double cosineHarmonicCoefficient,
        const double sineHarmonicCoefficient,
        const double legendrePolynomial,
        const double legendrePolynomialDerivative,
        const double legendrePolynomialSecondDerivative,
        Eigen::Matrix3d& sphericalHessian );

void computePotentialSphericalHessian(
        const Eigen::Vector3d& sphericalPosition,
        const double referenceRadius,
        const double preMultiplier,
        const int degree,
        const int order,
        const double cosineHarmonicCoefficient,
        const double sineHarmonicCoefficient,
        const double legendrePolynomial,
        const double legendrePolynomialDerivative,
        const double legendrePolynomialSecondDerivative,
        Eigen::Matrix3d& sphericalHessian );

void computePotentialSphericalHessian(
        const Eigen::Vector3d& sphericalPosition,
        const double preMultiplier,
        const int degree,
        const int order,
        const double cosineHarmonicCoefficient,
        const double sineHarmonicCoefficient,
        const boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache,
        Eigen::Matrix3d& sphericalHessian );

Eigen::Matrix3d computeCumulativeSphericalHessian(
        const Eigen::Vector3d& sphericalPosition,
        const double referenceRadius,
        const double gravitionalParameter,
        const Eigen::MatrixXd cosineHarmonicCoefficients,
        const Eigen::MatrixXd sineHarmonicCoefficients,
        const boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > shCache );

Eigen::Matrix3d computePartialDerivativeOfBodyFixedSphericalHarmonicAcceleration(
        const Eigen::Vector3d& cartesianPosition,
        const Eigen::Vector3d& sphericalPosition,
        const double referenceRadius,
        const double gravitionalParameter,
        const Eigen::MatrixXd cosineHarmonicCoefficients,
        const Eigen::MatrixXd sineHarmonicCoefficients,
        const boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > shCache,
        const Eigen::Vector3d& sphericalPotentialGradient,
        const Eigen::Matrix3d& sphericalToCartesianGradientMatrix );

Eigen::Matrix3d computePartialDerivativeOfBodyFixedSphericalHarmonicAcceleration(
        const Eigen::Vector3d& cartesianPosition,
        const double referenceRadius,
        const double gravitionalParameter,
        const Eigen::MatrixXd cosineHarmonicCoefficients,
        const Eigen::MatrixXd sineHarmonicCoefficients,
        const boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > shCache );

}

}

}
#endif // SPHERICALHARMONICPARTIALFUNCTIONS_H
