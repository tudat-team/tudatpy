#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"

#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityModel.h"


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
//    std::cout<<"Input "<<degree<<" "<<order<<" "<<radiusPowerTerm<<" "<<preMultiplier<<" "<<
//               cosineHarmonicCoefficient<<" "<<
//                      sineHarmonicCoefficient<<" "<<
//                      legendrePolynomial<<" "<<
//                      legendrePolynomialDerivative<<" "<<
//                      legendrePolynomialSecondDerivative<<" "<<std::endl;
    sphericalHessian.setZero( );

    sphericalHessian( 0, 0 ) += static_cast< double >( degree + 1 ) * static_cast< double >( degree + 2 ) /
            ( distance * distance ) * legendrePolynomial *
            (  cosineHarmonicCoefficient * cosineOfOrderLongitude + sineHarmonicCoefficient * sineOfOrderLongitude );
    sphericalHessian( 1, 0 ) += -static_cast< double >( degree + 1 ) / distance * cosineOfLatitude *
            legendrePolynomialDerivative *
            (  cosineHarmonicCoefficient * cosineOfOrderLongitude + sineHarmonicCoefficient * sineOfOrderLongitude );
    sphericalHessian( 2, 0 ) += -static_cast< double >( order ) * static_cast< double >( degree + 1 ) / distance *
            legendrePolynomial *
            ( -cosineHarmonicCoefficient * sineOfOrderLongitude + sineHarmonicCoefficient * cosineOfOrderLongitude );

    sphericalHessian( 0, 1 ) = sphericalHessian( 1, 0 );
    sphericalHessian( 1, 1 ) = ( cosineOfLatitude * cosineOfLatitude * legendrePolynomialSecondDerivative -
                                 sineOfLatitude * legendrePolynomialDerivative ) *
           ( cosineHarmonicCoefficient * cosineOfOrderLongitude + sineHarmonicCoefficient * sineOfOrderLongitude );
    sphericalHessian( 2, 1 ) = static_cast< double >( order ) * cosineOfLatitude * legendrePolynomialDerivative *
            ( -cosineHarmonicCoefficient * sineOfOrderLongitude + sineHarmonicCoefficient * cosineOfOrderLongitude );

    sphericalHessian( 0, 2 ) = sphericalHessian( 2, 0 );
    sphericalHessian( 1, 2 ) = sphericalHessian( 2, 1 );
    sphericalHessian( 2, 2 ) += static_cast< double >( order ) * static_cast< double >( order ) * legendrePolynomial *
            (  -cosineHarmonicCoefficient * cosineOfOrderLongitude - sineHarmonicCoefficient * sineOfOrderLongitude );


    sphericalHessian *= preMultiplier * radiusPowerTerm;

//    std::cout<<"Output:  "<<std::endl<<sphericalHessian<<std::endl<<std::endl;
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
    std::cout<<"Deg: "<<degree + 1<<" "<<shCache->getReferenceRadiusRatioPowers( degree + 1 )<<std::endl;
    computePotentialSphericalHessian(
                sphericalPosition( 0 )  , shCache->getReferenceRadiusRatioPowers( degree + 1 ),
                shCache->getCosineOfMultipleLongitude( order ), shCache->getSineOfMultipleLongitude( order ),
                shCache->getLegendreCache( )->getCurrentPolynomialParameterComplement( ),
                shCache->getLegendreCache( )->getCurrentPolynomialParameter( ),
                preMultiplier,
                degree, order, cosineHarmonicCoefficient, sineHarmonicCoefficient,
                shCache->getLegendreCache( )->getLegendrePolynomial( degree, order ),
                shCache->getLegendreCache( )->getLegendrePolynomialDerivative( degree, order ),
                shCache->getLegendreCache( )->getLegendrePolynomialSecondDerivative( degree, order ),
                sphericalHessian );
}

Eigen::Matrix3d computePartialDerivativeOfBodyFixedSphericalHarmonicAcceleration(
        const Eigen::Vector3d& cartesianPosition,
        const Eigen::Vector3d& sphericalPosition,
        const double referenceRadius,
        const double gravitionalParameter,
        const Eigen::MatrixXd cosineHarmonicCoefficients,
        const Eigen::MatrixXd sineHarmonicCoefficients,
        const boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > shCache,
        const Eigen::Vector3d& sphericalPotentialGradient,
        const Eigen::Matrix3d& sphericalToCartesianGradientMatrix )
{
    double preMultiplier = gravitionalParameter / referenceRadius;

    Eigen::Matrix3d sphericalHessian, sphericalHessianTerm;

    sphericalHessian.setZero( );
    for( unsigned int i = 0; i < cosineHarmonicCoefficients.rows( ); i++ )
    {
        for( unsigned int j = 0; ( j <= i && j < cosineHarmonicCoefficients.cols( ) ); j++ )
        {
            computePotentialSphericalHessian(
                        sphericalPosition, preMultiplier, i, j, cosineHarmonicCoefficients( i, j ),
                        sineHarmonicCoefficients( i, j ), shCache, sphericalHessianTerm );
            sphericalHessian += sphericalHessianTerm;
        }
    }

    Eigen::Matrix3d accelerationPartial =
            sphericalToCartesianGradientMatrix * sphericalHessian * sphericalToCartesianGradientMatrix.inverse( );

    accelerationPartial += coordinate_conversions::getDerivativeOfSphericalToCartesianGradient(
                cartesianPosition, sphericalPotentialGradient );

    std::cout<<"Contributions: "<<std::endl<<
               sphericalToCartesianGradientMatrix<<std::endl<<std::endl<<
               sphericalToCartesianGradientMatrix.inverse( )<<std::endl<<std::endl<<
               sphericalHessian<<std::endl<<std::endl<<
               sphericalToCartesianGradientMatrix * sphericalHessian * sphericalToCartesianGradientMatrix.inverse( )<<std::endl<<std::endl<<
               coordinate_conversions::getDerivativeOfSphericalToCartesianGradient(
                               cartesianPosition, sphericalPotentialGradient )<<std::endl;



    return accelerationPartial;
}

Eigen::Matrix3d computePartialDerivativeOfBodyFixedSphericalHarmonicAcceleration(
        const Eigen::Vector3d& cartesianPosition,
        const double referenceRadius,
        const double gravitionalParameter,
        const Eigen::MatrixXd cosineHarmonicCoefficients,
        const Eigen::MatrixXd sineHarmonicCoefficients,
        const boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > shCache )
{

    Eigen::Matrix3d gradientTransformationMatrix =
            coordinate_conversions::getSphericalToCartesianGradientMatrix( cartesianPosition );

    Eigen::Vector3d sphericalPosition =
            coordinate_conversions::convertCartesianToSpherical( cartesianPosition );
    sphericalPosition( 1 ) = mathematical_constants::PI / 2.0 - sphericalPosition( 1 );

    Eigen::Vector3d sphericalPotentialGradient = gradientTransformationMatrix.transpose( ) *
            gravitation::computeGeodesyNormalizedGravitationalAccelerationSum(
                cartesianPosition, gravitionalParameter, referenceRadius, cosineHarmonicCoefficients, sineHarmonicCoefficients,
                shCache );

    std::cout<<"Spherical: "<<sphericalPosition<<std::endl<<std::endl;
    std::cout<<"Inverse: "<<std::endl<<
               gradientTransformationMatrix.inverse( )<<std::endl<<"Transpose: "<<std::endl<<gradientTransformationMatrix.transpose( )<<std::endl;

    return computePartialDerivativeOfBodyFixedSphericalHarmonicAcceleration(
                cartesianPosition, sphericalPosition, referenceRadius, gravitionalParameter, cosineHarmonicCoefficients,
                sineHarmonicCoefficients, shCache, sphericalPotentialGradient, gradientTransformationMatrix );
}

}

}

}
