/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"

#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityModel.h"


#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/sphericalHarmonicPartialFunctions.h"

namespace tudat
{

namespace acceleration_partials
{

using namespace orbital_element_conversions;

//! Function to compute the spherical Hessian of a single term of a spherical harmonic potential
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

//! Function to compute the spherical Hessian of a single term of a spherical harmonic potential
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

//! Function to compute the spherical Hessian of a single term of a spherical harmonic potential
void computePotentialSphericalHessian(
        const Eigen::Vector3d& sphericalPosition,
        const double preMultiplier,
        const int degree,
        const int order,
        const double cosineHarmonicCoefficient,
        const double sineHarmonicCoefficient,
        const std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache,
        Eigen::Matrix3d& sphericalHessian )
{
    computePotentialSphericalHessian(
                sphericalPosition( 0 )  , sphericalHarmonicsCache->getReferenceRadiusRatioPowers( degree + 1 ),
                sphericalHarmonicsCache->getCosineOfMultipleLongitude( order ),
                sphericalHarmonicsCache->getSineOfMultipleLongitude( order ),
                sphericalHarmonicsCache->getLegendreCache( )->getCurrentPolynomialParameterComplement( ),
                sphericalHarmonicsCache->getLegendreCache( )->getCurrentPolynomialParameter( ),
                preMultiplier,
                degree, order, cosineHarmonicCoefficient, sineHarmonicCoefficient,
                sphericalHarmonicsCache->getLegendreCache( )->getLegendrePolynomial( degree, order ),
                sphericalHarmonicsCache->getLegendreCache( )->getLegendrePolynomialDerivative( degree, order ),
                sphericalHarmonicsCache->getLegendreCache( )->getLegendrePolynomialSecondDerivative( degree, order ),
                sphericalHessian );
}

//! Function to compute the spherical Hessian of a full spherical harmonic potential
Eigen::Matrix3d computeCumulativeSphericalHessian(
        const Eigen::Vector3d& sphericalPosition,
        const double referenceRadius,
        const double gravitionalParameter,
        const Eigen::MatrixXd cosineHarmonicCoefficients,
        const Eigen::MatrixXd sineHarmonicCoefficients,
        const std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache )
{
    double preMultiplier = gravitionalParameter / referenceRadius;

    Eigen::Matrix3d sphericalHessian, sphericalHessianTerm;

    sphericalHessian.setZero( );
    for( int i = 0; i < cosineHarmonicCoefficients.rows( ); i++ )
    {
        for( int j = 0; ( j <= i && j < cosineHarmonicCoefficients.cols( ) ); j++ )
        {
            computePotentialSphericalHessian(
                        sphericalPosition, preMultiplier, i, j, cosineHarmonicCoefficients( i, j ),
                        sineHarmonicCoefficients( i, j ), sphericalHarmonicsCache, sphericalHessianTerm );
            sphericalHessian += sphericalHessianTerm;
        }
    }
    return sphericalHessian;

}

//! Calculate partial of spherical harmonic acceleration w.r.t. position of body undergoing acceleration
//! (in the body-fixed frame)
Eigen::Matrix3d computePartialDerivativeOfBodyFixedSphericalHarmonicAcceleration(
        const Eigen::Vector3d& cartesianPosition,
        const Eigen::Vector3d& sphericalPosition,
        const double referenceRadius,
        const double gravitionalParameter,
        const Eigen::MatrixXd cosineHarmonicCoefficients,
        const Eigen::MatrixXd sineHarmonicCoefficients,
        const std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache,
        const Eigen::Vector3d& sphericalPotentialGradient,
        const Eigen::Matrix3d& sphericalToCartesianGradientMatrix )
{
    // Compute Hessian in spherical coordinates.
    Eigen::Matrix3d sphericalHessian = computeCumulativeSphericalHessian(
                sphericalPosition, referenceRadius, gravitionalParameter, cosineHarmonicCoefficients,
                sineHarmonicCoefficients, sphericalHarmonicsCache );

    // Convert to Cartesian Hessian
    Eigen::Matrix3d accelerationPartial =
            sphericalToCartesianGradientMatrix * sphericalHessian * sphericalToCartesianGradientMatrix.transpose( );

    // Add effect of direct change in rotation matrix
    accelerationPartial += coordinate_conversions::getDerivativeOfSphericalToCartesianGradient(
                sphericalPotentialGradient, cartesianPosition );

    return accelerationPartial;
}

//! Calculate partial of spherical harmonic acceleration w.r.t. position of body undergoing acceleration
//! (in the body-fixed frame)
Eigen::Matrix3d computePartialDerivativeOfBodyFixedSphericalHarmonicAcceleration(
        const Eigen::Vector3d& cartesianPosition,
        const double referenceRadius,
        const double gravitionalParameter,
        const Eigen::MatrixXd cosineHarmonicCoefficients,
        const Eigen::MatrixXd sineHarmonicCoefficients,
        const std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache )
{
    // Compute spherical position.
    Eigen::Vector3d sphericalPosition =
            coordinate_conversions::convertCartesianToSpherical( cartesianPosition );
    sphericalPosition( 1 ) = mathematical_constants::PI / 2.0 - sphericalPosition( 1 );

    // Compute spherical to Cartesian gradient transformation.
    Eigen::Matrix3d gradientTransformationMatrix =
            coordinate_conversions::getSphericalToCartesianGradientMatrix( cartesianPosition );

    // Compute spherical gradient.
    std::map< std::pair< int, int >, Eigen::Vector3d > dummyMap;
    Eigen::Vector3d sphericalPotentialGradient = gradientTransformationMatrix.inverse( ) *
            gravitation::computeGeodesyNormalizedGravitationalAccelerationSum(
                cartesianPosition, gravitionalParameter, referenceRadius, cosineHarmonicCoefficients,
                sineHarmonicCoefficients, sphericalHarmonicsCache, dummyMap );

    return computePartialDerivativeOfBodyFixedSphericalHarmonicAcceleration(
                cartesianPosition, sphericalPosition, referenceRadius, gravitionalParameter, cosineHarmonicCoefficients,
                sineHarmonicCoefficients, sphericalHarmonicsCache, sphericalPotentialGradient,
                gradientTransformationMatrix );
}

//! Calculate partial of spherical harmonic acceleration w.r.t. a set of cosine coefficients
void calculateSphericalHarmonicGravityWrtCCoefficients(
        const Eigen::Vector3d& sphericalPosition,
        const double referenceRadius,
        const double gravitionalParameter,
        const std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache,
        const std::vector< std::pair< int, int > >& blockIndices,
        const Eigen::Matrix3d& sphericalToCartesianGradientMatrix,
        const Eigen::Matrix3d& bodyFixedToIntegrationFrame,
        Eigen::MatrixXd& partialsMatrix )
{
    double preMultiplier = gravitionalParameter / referenceRadius;
    const std::shared_ptr< basic_mathematics::LegendreCache > legendreCache = sphericalHarmonicsCache->getLegendreCache( );

    int degree, order;
    for( unsigned int i = 0; i < blockIndices.size( ); i++ )
    {
        degree = blockIndices.at( i ).first;
        order = blockIndices.at( i ).second;

        // Calculate and set partial of current degree and order.
        partialsMatrix.block( 0, i, 3, 1 ) =
                basic_mathematics::computePotentialGradient(
                    sphericalPosition( radiusIndex ),
                    sphericalHarmonicsCache->getReferenceRadiusRatioPowers( degree + 1 ),
                    sphericalHarmonicsCache->getCosineOfMultipleLongitude( order ),
                    sphericalHarmonicsCache->getSineOfMultipleLongitude( order ),
                    sphericalHarmonicsCache->getLegendreCache( )->getCurrentPolynomialParameterComplement( ),
                    preMultiplier, degree, order,
                    1.0, 0.0, legendreCache->getLegendrePolynomial( degree, order ),
                    legendreCache->getLegendrePolynomialDerivative( degree, order ) );

    }

    // Transform partials to Cartesian position and integration frame.
    partialsMatrix = bodyFixedToIntegrationFrame * sphericalToCartesianGradientMatrix * partialsMatrix;
}

//! Calculate partial of spherical harmonic acceleration w.r.t. a set of sine coefficients
void calculateSphericalHarmonicGravityWrtSCoefficients(
        const Eigen::Vector3d& sphericalPosition,
        const double referenceRadius,
        const double gravitionalParameter,
        const std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache,
        const std::vector< std::pair< int, int > >& blockIndices,
        const Eigen::Matrix3d& sphericalToCartesianGradientMatrix,
        const Eigen::Matrix3d& bodyFixedToIntegrationFrame,
        Eigen::MatrixXd& partialsMatrix )
{
    double preMultiplier = gravitionalParameter / referenceRadius;
    const std::shared_ptr< basic_mathematics::LegendreCache > legendreCache = sphericalHarmonicsCache->getLegendreCache( );

    int degree, order;
    for( unsigned int i = 0; i < blockIndices.size( ); i++ )
    {
        degree = blockIndices.at( i ).first;
        order = blockIndices.at( i ).second;

        // Calculate and set partial of current degree and order.
        partialsMatrix.block( 0, i, 3, 1 ) =
                basic_mathematics::computePotentialGradient(
                    sphericalPosition( radiusIndex ),
                    sphericalHarmonicsCache->getReferenceRadiusRatioPowers( degree + 1 ),
                    sphericalHarmonicsCache->getCosineOfMultipleLongitude( order ),
                    sphericalHarmonicsCache->getSineOfMultipleLongitude( order ),
                    sphericalHarmonicsCache->getLegendreCache( )->getCurrentPolynomialParameterComplement( ),
                    preMultiplier, degree, order,
                    0.0, 1.0, legendreCache->getLegendrePolynomial( degree, order ),
                    legendreCache->getLegendrePolynomialDerivative( degree, order ) );

    }

    // Transform partials to Cartesian position and integration frame.
    partialsMatrix = bodyFixedToIntegrationFrame * sphericalToCartesianGradientMatrix * partialsMatrix;
}

}

}
