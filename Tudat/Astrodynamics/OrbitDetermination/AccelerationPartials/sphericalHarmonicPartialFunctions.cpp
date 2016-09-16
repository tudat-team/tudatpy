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
        const boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache,
        Eigen::Matrix3d& sphericalHessian )
{
    computePotentialSphericalHessian(
                sphericalPosition( 0 )  , sphericalHarmonicsCache->getReferenceRadiusRatioPowers( degree + 1 ),
                sphericalHarmonicsCache->getCosineOfMultipleLongitude( order ), sphericalHarmonicsCache->getSineOfMultipleLongitude( order ),
                sphericalHarmonicsCache->getLegendreCache( )->getCurrentPolynomialParameterComplement( ),
                sphericalHarmonicsCache->getLegendreCache( )->getCurrentPolynomialParameter( ),
                preMultiplier,
                degree, order, cosineHarmonicCoefficient, sineHarmonicCoefficient,
                sphericalHarmonicsCache->getLegendreCache( )->getLegendrePolynomial( degree, order ),
                sphericalHarmonicsCache->getLegendreCache( )->getLegendrePolynomialDerivative( degree, order ),
                sphericalHarmonicsCache->getLegendreCache( )->getLegendrePolynomialSecondDerivative( degree, order ),
                sphericalHessian );
}

Eigen::Matrix3d computeCumulativeSphericalHessian(
        const Eigen::Vector3d& sphericalPosition,
        const double referenceRadius,
        const double gravitionalParameter,
        const Eigen::MatrixXd cosineHarmonicCoefficients,
        const Eigen::MatrixXd sineHarmonicCoefficients,
        const boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache )
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
                        sineHarmonicCoefficients( i, j ), sphericalHarmonicsCache, sphericalHessianTerm );
            sphericalHessian += sphericalHessianTerm;
        }
    }
    return sphericalHessian;

}

Eigen::Matrix3d computePartialDerivativeOfBodyFixedSphericalHarmonicAcceleration(
        const Eigen::Vector3d& cartesianPosition,
        const Eigen::Vector3d& sphericalPosition,
        const double referenceRadius,
        const double gravitionalParameter,
        const Eigen::MatrixXd cosineHarmonicCoefficients,
        const Eigen::MatrixXd sineHarmonicCoefficients,
        const boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache,
        const Eigen::Vector3d& sphericalPotentialGradient,
        const Eigen::Matrix3d& sphericalToCartesianGradientMatrix )
{
    Eigen::Matrix3d sphericalHessian = computeCumulativeSphericalHessian(
                sphericalPosition, referenceRadius, gravitionalParameter, cosineHarmonicCoefficients,
                sineHarmonicCoefficients, sphericalHarmonicsCache );

    Eigen::Matrix3d accelerationPartial =
            sphericalToCartesianGradientMatrix * sphericalHessian * sphericalToCartesianGradientMatrix.transpose( );

    accelerationPartial += coordinate_conversions::getDerivativeOfSphericalToCartesianGradient(
                sphericalPotentialGradient, cartesianPosition );
    return accelerationPartial;
}

Eigen::Matrix3d computePartialDerivativeOfBodyFixedSphericalHarmonicAcceleration(
        const Eigen::Vector3d& cartesianPosition,
        const double referenceRadius,
        const double gravitionalParameter,
        const Eigen::MatrixXd cosineHarmonicCoefficients,
        const Eigen::MatrixXd sineHarmonicCoefficients,
        const boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache )
{

    Eigen::Matrix3d gradientTransformationMatrix =
            coordinate_conversions::getSphericalToCartesianGradientMatrix( cartesianPosition );

    Eigen::Vector3d sphericalPosition =
            coordinate_conversions::convertCartesianToSpherical( cartesianPosition );
    sphericalPosition( 1 ) = mathematical_constants::PI / 2.0 - sphericalPosition( 1 );

    Eigen::Vector3d sphericalPotentialGradient = gradientTransformationMatrix.inverse( ) *
            gravitation::computeGeodesyNormalizedGravitationalAccelerationSum(
                cartesianPosition, gravitionalParameter, referenceRadius, cosineHarmonicCoefficients, sineHarmonicCoefficients,
                sphericalHarmonicsCache );

    return computePartialDerivativeOfBodyFixedSphericalHarmonicAcceleration(
                cartesianPosition, sphericalPosition, referenceRadius, gravitionalParameter, cosineHarmonicCoefficients,
                sineHarmonicCoefficients, sphericalHarmonicsCache, sphericalPotentialGradient, gradientTransformationMatrix );
}

void calculateSphericalHarmonicGravityWrtCCoefficients(
        const Eigen::Vector3d& sphericalPosition,
        const double referenceRadius,
        const double gravitionalParameter,
        const boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache,
        const std::map< int, std::pair< int, int > >& blockIndices,
        const Eigen::Matrix3d& sphericalToCartesianGradientMatrix,
        const Eigen::Matrix3d& bodyFixedToIntegrationFrame,
        Eigen::MatrixXd& partialsMatrix )
{
    double preMultiplier = gravitionalParameter / referenceRadius;
    const boost::shared_ptr< basic_mathematics::LegendreCache > legendreCache = sphericalHarmonicsCache->getLegendreCache( );

    int currentIndex = 0;
    int degree;
    for( std::map< int, std::pair< int, int > >::const_iterator blockIterator = blockIndices.begin( );
         blockIterator != blockIndices.end( ); blockIterator++ )
    {
        degree = blockIterator->first;

        // Iterate over all required orders in current degree.
        for( int order = blockIterator->second.first; order < blockIterator->second.second + blockIterator->second.first; order++ )
        {

            // Calculate and set partial of current degree and order.
            partialsMatrix.block( 0, currentIndex, 3, 1 ) =
                    basic_mathematics::computePotentialGradient(
                        sphericalPosition( radiusIndex ),
                        sphericalHarmonicsCache->getReferenceRadiusRatioPowers( degree + 1 ),
                        sphericalHarmonicsCache->getCosineOfMultipleLongitude( order ),
                        sphericalHarmonicsCache->getSineOfMultipleLongitude( order ),
                        sphericalHarmonicsCache->getLegendreCache( )->getCurrentPolynomialParameterComplement( ),
                        preMultiplier, degree, order,
                        1.0, 0.0, legendreCache->getLegendrePolynomial( degree, order ),
                        legendreCache->getLegendrePolynomialDerivative( degree, order ) );

            currentIndex++;
        }
    }

    partialsMatrix = bodyFixedToIntegrationFrame * sphericalToCartesianGradientMatrix * partialsMatrix;
}

void calculateSphericalHarmonicGravityWrtSCoefficients(
        const Eigen::Vector3d& sphericalPosition,
        const double referenceRadius,
        const double gravitionalParameter,
        const boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache,
        const std::map< int, std::pair< int, int > >& blockIndices,
        const Eigen::Matrix3d& sphericalToCartesianGradientMatrix,
        const Eigen::Matrix3d& bodyFixedToIntegrationFrame,
        Eigen::MatrixXd& partialsMatrix )
{
    double preMultiplier = gravitionalParameter / referenceRadius;
    const boost::shared_ptr< basic_mathematics::LegendreCache > legendreCache = sphericalHarmonicsCache->getLegendreCache( );

    int currentIndex = 0;
    int degree;
    for( std::map< int, std::pair< int, int > >::const_iterator blockIterator = blockIndices.begin( );
         blockIterator != blockIndices.end( ); blockIterator++ )
    {
        degree = blockIterator->first;

        // Iterate over all required orders in current degree.
        for( int order = blockIterator->second.first; order < blockIterator->second.second + blockIterator->second.first; order++ )
        {

            // Calculate and set partial of current degree and order.
            partialsMatrix.block( 0, currentIndex, 3, 1 ) =
                    basic_mathematics::computePotentialGradient(
                        sphericalPosition( radiusIndex ),
                        sphericalHarmonicsCache->getReferenceRadiusRatioPowers( degree + 1 ),
                        sphericalHarmonicsCache->getCosineOfMultipleLongitude( order ),
                        sphericalHarmonicsCache->getSineOfMultipleLongitude( order ),
                        sphericalHarmonicsCache->getLegendreCache( )->getCurrentPolynomialParameterComplement( ),
                        preMultiplier, degree, order,
                        0.0, 1.0, legendreCache->getLegendrePolynomial( degree, order ),
                        legendreCache->getLegendrePolynomialDerivative( degree, order ) );

            currentIndex++;
        }
    }

    partialsMatrix = bodyFixedToIntegrationFrame * sphericalToCartesianGradientMatrix * partialsMatrix;
}

}

}

}
