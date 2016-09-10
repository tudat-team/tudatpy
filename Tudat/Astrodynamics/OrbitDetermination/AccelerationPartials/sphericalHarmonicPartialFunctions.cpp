#include "Tudat/Mathematics/BasicMathematics/sphericalHarmonics.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/sphericalHarmonicPartialFunctions.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{

using namespace orbital_element_conversions;

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
        Eigen::Matrix3d& sphericalHessian )
{
    sphericalHessian.setZero( );

    sphericalHessian( 0, 0 ) += static_cast< double >( degree + 1 ) * static_cast< double >( degree + 2 ) /
            ( distance * distance ) * legendrePolynomial *
            (  cosineHarmonicCoefficient * cosineOfOrderLongitude + sineHarmonicCoefficient * sineOfOrderLongitude );
    sphericalHessian( 1, 0 ) += -static_cast< double >( degree + 1 ) / distance * cosineOfLatitude *
            legendrePolynomialDerivative *
            (  cosineHarmonicCoefficient * cosineOfOrderLongitude + sineHarmonicCoefficient * sineOfOrderLongitude );
    sphericalHessian( 2, 0 ) += static_cast< double >( order ) * static_cast< double >( degree + 1 ) / distance *
            legendrePolynomial *
            ( -cosineHarmonicCoefficient * sineOfOrderLongitude + sineHarmonicCoefficient * cosineOfOrderLongitude );

    sphericalHessian( 0, 1 ) = sphericalHessian( 0, 1 );
    sphericalHessian( 1, 1 ) = ( cosineOfLatitude * legendrePolynomialSecondDerivative -
                                 sineOfLatitude * legendrePolynomialDerivative ) *
           ( cosineHarmonicCoefficient * cosineOfOrderLongitude + sineHarmonicCoefficient * sineOfOrderLongitude );
    sphericalHessian( 2, 1 ) = static_cast< double >( order ) * cosineOfLatitude * legendrePolynomialDerivative *
            ( -cosineHarmonicCoefficient * sineOfOrderLongitude + sineHarmonicCoefficient * cosineOfOrderLongitude );

    sphericalHessian( 0, 2 ) = sphericalHessian( 2, 0 );
    sphericalHessian( 1, 2 ) = sphericalHessian( 2, 1 );
    sphericalHessian( 2, 2 ) += static_cast< double >( order ) * static_cast< double >( order ) * legendrePolynomial *
            (  -cosineHarmonicCoefficient * cosineOfOrderLongitude - sineHarmonicCoefficient * sineOfOrderLongitude );

    sphericalHessian *= preMultiplier * radiusPowerTerm;
}

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
        Eigen::Matrix3d& sphericalHessian )
{
    computePotentialSphericalHessian(
                sphericalPosition( radiusIndex ),
                basic_mathematics::raiseToIntegerPower
                ( referenceRadius / sphericalPosition( radiusIndex ), static_cast< double >( degree ) + 1.0 ),
                std::cos( static_cast< double >( order ) * sphericalPosition( longitudeIndex ) ),
                std::sin( static_cast< double >( order ) * sphericalPosition( longitudeIndex ) ),
                std::cos( sphericalPosition( latitudeIndex ) ), std::sin( sphericalPosition( latitudeIndex ) ),
                preMultiplier, degree, order, cosineHarmonicCoefficient, sineHarmonicCoefficient,
                legendrePolynomial,legendrePolynomialDerivative, legendrePolynomialSecondDerivative, sphericalHessian );
}

void computePotentialSphericalHessian(
        const Eigen::Vector3d& sphericalPosition,
        const double preMultiplier,
        const int degree,
        const int order,
        const double cosineHarmonicCoefficient,
        const double sineHarmonicCoefficient,
        const boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > shCache,
        Eigen::Matrix3d& sphericalHessian )
{
    computePotentialSphericalHessian(
                sphericalPosition, shCache->getReferenceRadiusRatioPowers( degree + 1 ), preMultiplier,
                degree, order, cosineHarmonicCoefficient, sineHarmonicCoefficient,
                shCache->getLegendreCache( )->getLegendrePolynomial( degree, order ),
                shCache->getLegendreCache( )->getLegendrePolynomialDerivative( degree, order ),
                shCache->getLegendreCache( )->getLegendrePolynomialSecondDerivative( degree, order ),
                sphericalHessian );

}


}

}

}
