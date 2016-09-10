#ifndef SPHERICALHARMONICPARTIALFUNCTIONS_H
#define SPHERICALHARMONICPARTIALFUNCTIONS_H



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
        const boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > shCache,
        Eigen::Matrix3d& sphericalHessian );

}

}

}
#endif // SPHERICALHARMONICPARTIALFUNCTIONS_H
